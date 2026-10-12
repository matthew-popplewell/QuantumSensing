#=
Instrument description and focal-plane geometry.

Focal-plane convention: x (column) increases toward +RA (east on the sky) and
y (row) toward +Dec (north). Pixel centres are at integer coordinates, with
the boresight at the centre of the (sub)frame.
=#

"""
    Telescope(; aperture_diameter, focal_length, obstruction_ratio=0.0,
              throughput=1.0)

Telescope optics.

# Fields
- `aperture_diameter::Meters`: Primary mirror diameter.
- `focal_length::Meters`: Effective focal length.
- `obstruction_ratio::Float64`: Central obstruction, as a fraction of the
  primary **diameter**.
- `throughput::Float64`: Product of optical transmissions (mirrors,
  corrector, filter), 0-1.
"""
struct Telescope
    aperture_diameter::Meters
    focal_length::Meters
    obstruction_ratio::Float64
    throughput::Float64
end

function Telescope(;
    aperture_diameter,
    focal_length,
    obstruction_ratio = 0.0,
    throughput = 1.0
)
    return Telescope(
        to_unit(u"m", aperture_diameter),
        to_unit(u"m", focal_length),
        to_fraction(obstruction_ratio),
        to_fraction(throughput)
    )
end

"""
    collecting_area(telescope)

Unobstructed collecting area, `pi/4 D^2 (1 - obstruction_ratio^2)`.
"""
function collecting_area(t::Telescope)
    return (pi / 4) * t.aperture_diameter^2 * (1 - t.obstruction_ratio^2)
end

"""
    Detector(; pixel_pitch, npixels, quantum_efficiency, read_noise,
             dark_current, gain, full_well=Inf, bias=0.0, bit_depth=16)

Detector characteristics. Electron and ADU counts are plain numbers.

# Fields
- `pixel_pitch::Meters`: Pixel size.
- `npixels::Tuple{Int, Int}`: Full-frame size (rows, columns).
- `quantum_efficiency::Float64`: QE in the band, 0-1 (or `u"percent"`).
- `read_noise::Float64`: Read noise, e- RMS per pixel.
- `dark_current::PerSecond`: Dark current, e- per pixel per second.
- `gain::Float64`: e- per ADU.
- `full_well::Float64`: Full-well capacity, e- (saturation level).
- `bias::Float64`: Bias offset, ADU.
- `bit_depth::Int`: ADC bits.
"""
struct Detector
    pixel_pitch::Meters
    npixels::Tuple{Int, Int}
    quantum_efficiency::Float64
    read_noise::Float64
    dark_current::PerSecond
    gain::Float64
    full_well::Float64
    bias::Float64
    bit_depth::Int
end

function Detector(;
    pixel_pitch,
    npixels,
    quantum_efficiency,
    read_noise,
    dark_current,
    gain,
    full_well = Inf,
    bias = 0.0,
    bit_depth = 16
)
    return Detector(
        to_unit(u"m", pixel_pitch),
        Tuple{Int, Int}(npixels),
        to_fraction(quantum_efficiency),
        Float64(read_noise),
        to_unit(u"s^-1", dark_current),
        Float64(gain),
        Float64(full_well),
        Float64(bias),
        Int(bit_depth)
    )
end

"""
    Bandpass(; wavelength, width, zero_point, sun_magnitude)

Photometric band.

# Fields
- `wavelength::Meters`: Pivot wavelength (used for F_nu -> F_lambda and the
  mean photon energy).
- `width::Meters`: Effective width.
- `zero_point::SpectralFlux`: Flux density of a magnitude-0 source. Accepts
  F_lambda (e.g. `u"erg/s/cm^2/Å"`) or F_nu (e.g. `u"W/m^2/Hz"`); F_nu is
  converted with the pivot wavelength.
- `sun_magnitude::Float64`: Apparent magnitude of the Sun at 1 AU.
"""
struct Bandpass
    wavelength::Meters
    width::Meters
    zero_point::SpectralFlux
    sun_magnitude::Float64
end

function Bandpass(; wavelength, width, zero_point, sun_magnitude)
    lam = to_unit(u"m", wavelength)
    zp = if dimension(zero_point) == dimension(u"W/m^2/Hz")
        to_unit(u"W/m^3", zero_point * Unitful.c0 / lam^2)
    else
        to_unit(u"W/m^3", zero_point)
    end
    return Bandpass(lam, to_unit(u"m", width), zp, Float64(sun_magnitude))
end

"""
    Sensor(telescope, detector, band)

A telescope + detector + band combination.
"""
struct Sensor
    telescope::Telescope
    detector::Detector
    band::Bandpass
end

Sensor(; telescope, detector, band) = Sensor(telescope, detector, band)

"""
    pixel_scale(sensor)

Angular size of one pixel (arcsec).
"""
function pixel_scale(s::Sensor)
    ratio = uconvert(NoUnits, s.detector.pixel_pitch / s.telescope.focal_length)
    return uconvert(u"arcsec", ratio * u"rad")
end

"""
    field_of_view(sensor, npixels=sensor.detector.npixels)

Angular extent (rows, columns) of a frame of `npixels` (degrees).
"""
function field_of_view(s::Sensor, npixels::Tuple{Int, Int} = s.detector.npixels)
    return uconvert.(u"°", npixels .* pixel_scale(s))
end

"""
    Observation(; boresight, exposure_time, seeing_fwhm, tracking=:sidereal)

One exposure. Target and observer states passed alongside are taken at the
exposure start.

# Fields
- `boresight_ra::Radians`, `boresight_dec::Radians`: Pointing at exposure
  start (set with `boresight = (ra, dec)`, any angle units; see
  `point_at`).
- `exposure_time::Seconds`: Integration time.
- `seeing_fwhm::Radians`: Gaussian PSF FWHM.
- `tracking::Symbol`: Mount mode:
    - `:sidereal`: fixed on the stars (GEO objects streak at ~15"/s);
    - `:earth_fixed`: mount stopped (stars streak, GEO is ~stationary);
    - `:rate_track`: follows the target (point source).
"""
struct Observation
    boresight_ra::Radians
    boresight_dec::Radians
    exposure_time::Seconds
    seeing_fwhm::Radians
    tracking::Symbol
end

function Observation(;
    boresight,
    exposure_time,
    seeing_fwhm,
    tracking::Symbol = :sidereal
)
    tracking in (:sidereal, :earth_fixed, :rate_track) || throw(ArgumentError(
        "Unknown tracking mode $tracking. Use :sidereal, :earth_fixed " *
        "or :rate_track."
    ))
    return Observation(
        to_unit(u"rad", boresight[1]),
        to_unit(u"rad", boresight[2]),
        to_unit(u"s", exposure_time),
        to_unit(u"rad", seeing_fwhm),
        tracking
    )
end

"""
    point_at(target, observer)

Boresight `(ra, dec)` that centres `target` as seen from `observer`.
"""
point_at(target, observer) = radec(target, observer)

"""
    boresight_at(obs, target, observer, t)

Boresight `(ra, dec)` at time `t` after exposure start for the observation's
tracking mode.
"""
function boresight_at(obs::Observation, target, observer, t)
    t_s = to_unit(u"s", t)
    ra0, dec0 = obs.boresight_ra, obs.boresight_dec
    if obs.tracking == :sidereal
        return (ra0, dec0)
    elseif obs.tracking == :earth_fixed
        omega = uconvert(u"rad/s", EARTH_ROTATION_RATE)
        return (ra0 + omega * t_s, dec0)
    else  # :rate_track keeps the initial target offset from the boresight
        ra_a, dec_a = radec(target, observer)
        ra_b, dec_b = radec(drift(target, t_s), drift(observer, t_s))
        dra = rem2pi(ustrip(u"rad", ra_b - ra_a), RoundNearest) * u"rad"
        return (ra0 + dra, dec0 + (dec_b - dec_a))
    end
end

"""
    gnomonic(ra, dec, ra0, dec0)

Tangent-plane coordinates `(xi, eta)` (rad, plain Float64) of `(ra, dec)`
about `(ra0, dec0)`. Returns `(NaN, NaN)` for points >= 90 deg away.
"""
function gnomonic(ra, dec, ra0, dec0)
    cd, sd = cos(dec), sin(dec)
    cd0, sd0 = cos(dec0), sin(dec0)
    cda, sda = cos(ra - ra0), sin(ra - ra0)
    denom = sd0 * sd + cd0 * cd * cda
    denom > 0 || return (NaN, NaN)
    return (cd * sda / denom, (cd0 * sd - sd0 * cd * cda) / denom)
end

"""
    sky_pixel(sensor, ra, dec, boresight; frame_size=sensor.detector.npixels)

Continuous pixel coordinates `(x, y)` of a fixed sky position for a given
boresight `(ra0, dec0)`.
"""
function sky_pixel(
    s::Sensor,
    ra,
    dec,
    boresight;
    frame_size::Tuple{Int, Int} = s.detector.npixels
)
    xi, eta = gnomonic(ra, dec, boresight[1], boresight[2])
    plate = uconvert(NoUnits, s.telescope.focal_length / s.detector.pixel_pitch)
    nrows, ncols = frame_size
    return ((ncols + 1) / 2 + plate * xi, (nrows + 1) / 2 + plate * eta)
end

"""
    target_pixel(sensor, obs, target, observer, t=0u"s";
                 frame_size=sensor.detector.npixels)

Continuous pixel coordinates `(x, y)` = (column, row) of the target at time
`t` after exposure start, with the boresight at the frame centre.
"""
function target_pixel(
    s::Sensor,
    obs::Observation,
    target,
    observer,
    t = 0.0u"s";
    frame_size::Tuple{Int, Int} = s.detector.npixels
)
    t_s = to_unit(u"s", t)
    ra, dec = radec(drift(target, t_s), drift(observer, t_s))
    return sky_pixel(s, ra, dec, boresight_at(obs, target, observer, t_s);
                     frame_size = frame_size)
end

"""
    pixel_velocity(sensor, obs, target, observer; dt=0.5u"s")

Apparent target motion on the focal plane at exposure start, `(vx, vy)` in
pixels per second (plain Float64), by central difference of
`target_pixel`.
"""
function pixel_velocity(s::Sensor, obs::Observation, target, observer; dt = 0.5u"s")
    h = to_unit(u"s", dt)
    p_plus = target_pixel(s, obs, target, observer, h)
    p_minus = target_pixel(s, obs, target, observer, -h)
    scale = 2 * ustrip(u"s", h)
    return ((p_plus[1] - p_minus[1]) / scale, (p_plus[2] - p_minus[2]) / scale)
end
