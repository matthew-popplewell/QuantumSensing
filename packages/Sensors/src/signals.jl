#=
Radiometry and image formation.
=#

"""
    photon_energy(band)

Mean photon energy at the band's pivot wavelength.
"""
photon_energy(b::Bandpass) = uconvert(u"J", Unitful.h * Unitful.c0 / b.wavelength)

"""
    zero_point_photon_flux(band)

In-band photon flux of a magnitude-0 source (photons s^-1 m^-2).
"""
function zero_point_photon_flux(b::Bandpass)
    return uconvert(u"s^-1*m^-2", b.zero_point * b.width / photon_energy(b))
end

"""
    apparent_magnitude(band, target, sun, observer, radius, albedo)

In-band apparent magnitude of a Lambertian sphere. `sun` and `observer` are
states in the same frame/origin as `target`.
"""
function apparent_magnitude(b::Bandpass, target, sun, observer, radius, albedo)
    reflect = lambertian_extended_source(target, sun, observer, radius, albedo)
    d_sun = norm(CartesianState(sun).r - CartesianState(target).r)
    d_obs = norm(CartesianState(observer).r - CartesianState(target).r)
    flux_ratio = uconvert(NoUnits, reflect * (AU / d_sun)^2 / d_obs^2)
    return b.sun_magnitude - 2.5 * log10(flux_ratio)
end

"""
    electron_rate(sensor, magnitude)

Detected photo-electron rate (s^-1) from a point source of `magnitude`.
"""
function electron_rate(s::Sensor, magnitude::Real)
    return uconvert(
        u"s^-1",
        zero_point_photon_flux(s.band) * 10^(-0.4 * magnitude) *
        collecting_area(s.telescope) * s.telescope.throughput *
        s.detector.quantum_efficiency
    )
end

"""
    sky_electron_rate(sensor, sky_brightness)

Sky photo-electron rate per pixel (s^-1) for a surface brightness in
mag/arcsec^2.
"""
function sky_electron_rate(s::Sensor, sky_brightness::Real)
    per_arcsec2 = electron_rate(s, sky_brightness) / (1.0u"arcsec")^2
    return uconvert(u"s^-1", per_arcsec2 * pixel_scale(s)^2)
end

"""
    _pixel_weights(center, sigma, lo, hi)

Fraction of a 1-D Gaussian (in pixels) falling in each pixel `lo:hi`.
"""
function _pixel_weights(center, sigma, lo, hi)
    return [
        cdf(Normal(center, sigma), k + 0.5) - cdf(Normal(center, sigma), k - 0.5)
        for k in lo:hi
    ]
end

"""
    _render_track!(expected, track, electrons, sigma, nsteps)

Add `electrons` spread along `track(f)` (pixel position at exposure fraction
`f` in [0, 1]) with a Gaussian PSF of `sigma` px, sampled at `nsteps`
sub-steps. Light falling outside the frame is lost.
"""
function _render_track!(expected, track, electrons, sigma, nsteps)
    nrows, ncols = size(expected)
    box = ceil(Int, 5 * sigma)
    e_step = electrons / nsteps
    for k in 0:(nsteps - 1)
        xc, yc = track((k + 0.5) / nsteps)
        (isfinite(xc) && isfinite(yc)) || continue
        r_lo, r_hi = max(1, floor(Int, yc) - box), min(nrows, ceil(Int, yc) + box)
        c_lo, c_hi = max(1, floor(Int, xc) - box), min(ncols, ceil(Int, xc) + box)
        (r_lo > r_hi || c_lo > c_hi) && continue
        wy = _pixel_weights(yc, sigma, r_lo, r_hi)
        wx = _pixel_weights(xc, sigma, c_lo, c_hi)
        @views expected[r_lo:r_hi, c_lo:c_hi] .+= e_step .* (wy * wx')
    end
    return expected
end

# Sub-steps needed for a track: `per_pixel` samples per pixel of motion.
function _substeps(track, per_pixel)
    p0, p1 = track(0.0), track(1.0)
    len = hypot(p1[1] - p0[1], p1[2] - p0[2])
    return max(10, ceil(Int, per_pixel * (isfinite(len) ? len : 0.0)))
end

_unit_vector(ra, dec) = SVector(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec))

"""
    simulate_image(sensor, obs, site, et, target; radius, albedo,
                   sky=SkyConditions(), stars=true, star_limit=nothing,
                   frame_size=sensor.detector.npixels, substeps_per_pixel=4,
                   rng=Random.default_rng())

Simulate one exposure of a Lambertian-sphere target seen from a ground site.

The expected-electron image is built from:
- the target (penumbra-scaled, extinguished), smeared by the tracking mode
  and spread by a Gaussian seeing PSF;
- synthetic stars (Gaia DR3-calibrated counts), rendered the same way;
- the sky background (`sky_background`): airglow, light pollution, twilight,
  moonlight, zodiacal light and integrated starlight from stars fainter than
  `star_limit`;
- dark current.
A single Poisson draw is made per pixel, then Gaussian read noise,
full-well clipping, gain, bias and ADC clipping are applied (rounding to
integer ADU adds the g^2/12 quantisation noise).

The sky is evaluated at the boresight and taken as uniform over the frame.
Not modelled: flat field, blooming, cosmic rays, light-time and aberration.

# Arguments
- `sensor`, `obs`: Instrument and exposure.
- `site::GroundSite`, `et`: Observer and exposure start (ET seconds or a
  time string).
- `target`: Earth-centred J2000 target state at `et` (any units).
- `radius`, `albedo`: Lambertian-sphere size and reflectance.
- `sky::SkyConditions`: Environment settings.
- `stars::Bool`: Draw synthetic stars.
- `star_limit`: Faintest V magnitude drawn as individual stars. Fainter
  stars go into the diffuse starlight. Default: the magnitude giving
  `read_noise` electrons in the exposure, capped at 21.
- `frame_size`, `substeps_per_pixel`, `rng`: Read-out window, temporal
  sampling, random generator.

# Returns
A NamedTuple with `adu`, `expected_electrons`, `magnitude` (above the
atmosphere), `phase_angle`, `illumination`, `target_electrons`,
`start_pixel`, `end_pixel`, `background` (from `sky_background`),
`background_electrons` (per component, e-/px), `stars` (vector of
`(star, start_pixel, end_pixel, electrons)`) and `star_limit`.
"""
function simulate_image(
    s::Sensor,
    obs::Observation,
    site::GroundSite,
    et,
    target;
    radius,
    albedo,
    sky::SkyConditions = SkyConditions(),
    stars::Bool = true,
    star_limit = nothing,
    frame_size::Tuple{Int, Int} = s.detector.npixels,
    substeps_per_pixel::Integer = 4,
    rng::AbstractRNG = Random.default_rng()
)
    et = ephemeris_time(et)
    target = CartesianState(target)
    observer = site_state(site, et)
    sun = spice_state("SUN", et)
    det = s.detector
    nrows, ncols = frame_size
    t_exp = obs.exposure_time
    sigma = uconvert(NoUnits, obs.seeing_fwhm / pixel_scale(s)) /
            (2 * sqrt(2 * log(2)))

    # Pointing and atmosphere along the boresight.
    boresight = _unit_vector(obs.boresight_ra, obs.boresight_dec)
    zenith = horizon_coordinates(site, et, boresight).zenith
    transmission = atmospheric_transmission(sky, zenith)
    electrons_for(v_mag) = ustrip(NoUnits, electron_rate(s, v_mag) * t_exp)

    # Faintest star drawn individually.
    v_limit = star_limit === nothing ?
              min(21.0, 2.5 * log10(electrons_for(0.0) * transmission / det.read_noise)) :
              Float64(star_limit)

    expected = zeros(Float64, nrows, ncols)

    # Target.
    mag = apparent_magnitude(s.band, target, sun, observer, radius, albedo)
    illum = illumination_fraction(target, sun)
    target_e = electrons_for(mag) * illum * transmission
    target_track(f) = target_pixel(s, obs, target, observer, f * t_exp;
                                   frame_size = frame_size)
    _render_track!(expected, target_track, target_e, sigma,
                   _substeps(target_track, substeps_per_pixel))

    # Stars: draw over a cone covering the frame for the whole exposure.
    counts = star_counts(boresight)
    star_list = []
    if stars
        b0 = boresight_at(obs, target, observer, 0.0u"s")
        b1 = boresight_at(obs, target, observer, t_exp)
        u0, u1 = _unit_vector(b0...), _unit_vector(b1...)
        sweep = acos(clamp(dot(u0, u1), -1.0, 1.0)) * u"rad"
        half_diag = uconvert(u"rad", hypot(nrows, ncols) / 2 * pixel_scale(s))
        margin = uconvert(u"rad", 5 * sigma * pixel_scale(s))
        cone = half_diag + sweep / 2 + margin
        field, _ = synthetic_star_field(normalize(u0 + u1), cone, v_limit; rng = rng)
        for star in field
            track(f) = sky_pixel(s, star.ra, star.dec,
                                 boresight_at(obs, target, observer, f * t_exp);
                                 frame_size = frame_size)
            e = electrons_for(star.v_mag) * transmission
            _render_track!(expected, track, e, sigma, _substeps(track, substeps_per_pixel))
            push!(star_list, (star = star, start_pixel = track(0.0),
                              end_pixel = track(1.0), electrons = e))
        end
    end

    # Sky background and dark current.
    isl = integrated_starlight_s10(counts, stars ? v_limit : -Inf)
    bg = sky_background(s.band, sky, site, et, boresight; isl_s10 = isl)
    component_names = (:airglow, :light_pollution, :twilight, :moonlight,
                       :zodiacal, :starlight)
    bg_e = NamedTuple{component_names}(Tuple(
        isfinite(bg[k]) ? ustrip(NoUnits, sky_electron_rate(s, bg[k]) * t_exp) : 0.0
        for k in component_names
    ))
    dark_e = ustrip(NoUnits, det.dark_current * t_exp)
    expected .+= sum(values(bg_e)) + dark_e

    # Noise, saturation and digitisation.
    adc_max = 2^det.bit_depth - 1
    adu = Matrix{Int}(undef, nrows, ncols)
    read = Normal(0.0, det.read_noise)
    for i in eachindex(expected)
        electrons = rand(rng, Poisson(expected[i])) + rand(rng, read)
        electrons = min(electrons, det.full_well)
        adu[i] = clamp(round(Int, electrons / det.gain + det.bias), 0, adc_max)
    end

    return (
        adu = adu,
        expected_electrons = expected,
        magnitude = mag,
        phase_angle = solar_phase_angle(target, sun, observer) * u"rad",
        illumination = illum,
        target_electrons = target_e,
        start_pixel = target_track(0.0),
        end_pixel = target_track(1.0),
        background = bg,
        background_electrons = merge(bg_e, (dark = dark_e,)),
        stars = star_list,
        star_limit = v_limit
    )
end
