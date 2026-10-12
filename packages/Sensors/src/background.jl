#=
Sky background model (V band), following Robinson (2023), "Light curve
simulation and shape inversion for human-made space objects", MS thesis,
Purdue University, Sec. 2.5.5 and App. B, with the corrections noted inline
(search "CORRECTION").

Every component is returned as a V-band surface brightness (mag/arcsec^2)
along the line of sight; `sky_electron_rate` turns that into e-/px/s.

Atmospheric components (airglow, light pollution, twilight, moonlight) grow
with path length (Robinson multiplies by the relative airmass). Components
from outside the atmosphere (zodiacal light, integrated starlight, and the
target and stars) are attenuated by extinction.
=#

# ---------------------------------------------------------------------------
# Small interpolation helpers (regular or irregular ascending grids, clamped).
# ---------------------------------------------------------------------------
function _bracket(grid::AbstractVector, x)
    xc = clamp(x, first(grid), last(grid))
    i = clamp(searchsortedlast(grid, xc), 1, length(grid) - 1)
    w = (xc - grid[i]) / (grid[i + 1] - grid[i])
    return i, w
end

function _interp1(grid, values, x)
    i, w = _bracket(grid, x)
    return (1 - w) * values[i] + w * values[i + 1]
end

function _interp2(g1, g2, table, x1, x2)
    i, wi = _bracket(g1, x1)
    j, wj = _bracket(g2, x2)
    return (1 - wi) * (1 - wj) * table[i, j] + wi * (1 - wj) * table[i + 1, j] +
           (1 - wi) * wj * table[i, j + 1] + wi * wj * table[i + 1, j + 1]
end

function _interp3(g1, g2, g3, table, x1, x2, x3)
    k, wk = _bracket(g3, x3)
    lo = _interp2(g1, g2, view(table, :, :, k), x1, x2)
    hi = _interp2(g1, g2, view(table, :, :, k + 1), x1, x2)
    return (1 - wk) * lo + wk * hi
end

# Surface brightness <-> linear (relative to a 0 mag/arcsec^2 source).
_mag_to_flux(mu) = isfinite(mu) ? 10.0^(-0.4 * mu) : 0.0
_flux_to_mag(f) = f > 0 ? -2.5 * log10(f) : Inf

# S10: number of 10th-magnitude stars per square degree (Robinson Eq. 2.73).
s10_to_mag(b10) = b10 > 0 ? 27.78 - 2.5 * log10(b10) : Inf

# ---------------------------------------------------------------------------
# Atmosphere
# ---------------------------------------------------------------------------
"""
    relative_airmass(zenith_angle)

Pickering (2002) relative airmass, as used by Robinson Eq. 2.82. It stays
finite at the horizon (about 40), unlike sec(z).
"""
function relative_airmass(zenith_angle)
    h = 90.0 - ustrip(u"°", to_unit(u"°", zenith_angle))  # apparent altitude
    h = max(h, 0.0)
    return 1 / sind(h + 244 / (165 + 47 * h^1.1))
end

"""
    SkyConditions(; extinction=0.11, light_pollution=nothing,
                  airglow_scale=1.0, sodium_rayleigh=65.0,
                  mie_a1=75.0, mie_a2=0.015, mie_psi1=15u"°", mie_psi2=50u"°")

Environmental settings for the sky background.

# Fields
- `extinction::Float64`: V-band extinction coefficient, mag/airmass. The
  default is Robinson Table B.1 at 550 nm. CORRECTION: that table holds
  extinction coefficients k (mag/airmass), not transmissions. It is used as
  T = 10^(-0.4 k X).
- `light_pollution::Union{Nothing, Float64}`: Artificial zenith sky
  brightness at the site, V mag/arcsec^2 (e.g. from the Falchi et al. 2016
  World Atlas). `nothing` means no artificial light.
- `airglow_scale::Float64`: Multiplier on the Leinert (1998) airglow, e.g.
  for solar activity.
- `sodium_rayleigh::Float64`: Na D airglow intensity in rayleighs (Leinert
  Table 13: 30 in summer to 100 in winter).
- `mie_*`: Daniels (1977) Mie-scattering parameters (Robinson Eq. 2.88; the
  defaults are the centres of the recommended ranges).
"""
struct SkyConditions
    extinction::Float64
    light_pollution::Union{Nothing, Float64}
    airglow_scale::Float64
    sodium_rayleigh::Float64
    mie_a1::Float64
    mie_a2::Float64
    mie_psi1::Radians
    mie_psi2::Radians
end

function SkyConditions(;
    extinction = 0.11,
    light_pollution = nothing,
    airglow_scale = 1.0,
    sodium_rayleigh = 65.0,
    mie_a1 = 75.0,
    mie_a2 = 0.015,
    mie_psi1 = 15u"°",
    mie_psi2 = 50u"°"
)
    return SkyConditions(
        Float64(extinction),
        light_pollution === nothing ? nothing : Float64(light_pollution),
        Float64(airglow_scale),
        Float64(sodium_rayleigh),
        Float64(mie_a1),
        Float64(mie_a2),
        to_unit(u"rad", mie_psi1),
        to_unit(u"rad", mie_psi2)
    )
end

"""
    atmospheric_transmission(sky, zenith_angle)

Fraction of light from above the atmosphere that reaches the telescope,
`10^(-0.4 k X)`.
"""
function atmospheric_transmission(sky::SkyConditions, zenith_angle)
    return 10.0^(-0.4 * sky.extinction * relative_airmass(zenith_angle))
end

# ---------------------------------------------------------------------------
# Local horizon frame
# ---------------------------------------------------------------------------
"""
    horizon_coordinates(site, et, direction; frame="J2000")

Zenith angle and azimuth (from north, through east) of a unit `direction`
seen from `site`.
"""
function horizon_coordinates(
    site::GroundSite,
    et,
    direction::AbstractVector;
    frame::AbstractString = "J2000"
)
    up = site_zenith(site, et; frame = frame)
    pole = SVector{3}(pxform("IAU_EARTH", frame, ephemeris_time(et)) * [0.0, 0.0, 1.0])
    east = normalize(cross(pole, up))
    north = cross(up, east)
    u = normalize(SVector{3, Float64}(direction))
    zenith = acos(clamp(dot(u, up), -1.0, 1.0)) * u"rad"
    azimuth = mod2pi(atan(dot(u, east), dot(u, north))) * u"rad"
    return (zenith = zenith, azimuth = azimuth)
end

# ---------------------------------------------------------------------------
# Airglow (Leinert et al. 1998, Table 13, through the Bessell V curve)
# ---------------------------------------------------------------------------
# Bessell V transmission, SVO Filter Profile Service (Generic/Bessell.V):
# http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.V
const BESSELL_V_CURVE_ANGSTROM = collect(4700.0:100.0:7000.0)
const BESSELL_V_CURVE_T = [
    0.000, 0.030, 0.163, 0.458, 0.780, 0.967, 1.000, 0.973, 0.898, 0.792,
    0.684, 0.574, 0.461, 0.359, 0.270, 0.197, 0.135, 0.081, 0.045, 0.025,
    0.017, 0.013, 0.009, 0.000
]

_bessell_v(lambda_angstrom) =
    _interp1(BESSELL_V_CURVE_ANGSTROM, BESSELL_V_CURVE_T, lambda_angstrom)

"""
    airglow_zenith_photon_radiance(sky, band)

V-band zenith airglow photon radiance (photons s^-1 m^-2 sr^-1), from
Leinert et al. (1998) Table 13: [OI] 557.7 nm 250 R, Na D (`sodium_rayleigh`),
[OI] 630.0 nm 60 R, [OI] 636.4 nm 20 R, and a 0.3 R/Å pseudo-continuum.
Lines are weighted by the Bessell V transmission; the continuum by the band's
effective width. Robinson's airglow spectrum A(lambda) is not tabulated in
his thesis, so these values are used instead.
"""
function airglow_zenith_photon_radiance(sky::SkyConditions, band::Bandpass)
    lines = [(5577.0, 250.0), (5893.0, sky.sodium_rayleigh), (6300.0, 60.0),
             (6364.0, 20.0)]
    line_r = sum(r * _bessell_v(l) for (l, r) in lines)
    continuum_r = 0.3 * ustrip(u"Å", band.width)
    rayleigh = 1e10 / (4pi) * u"s^-1*m^-2*sr^-1"  # 1 R in photon radiance
    return sky.airglow_scale * (line_r + continuum_r) * rayleigh
end

"""
    photon_radiance_to_mag(band, radiance)

Convert an in-band photon radiance (photons s^-1 m^-2 sr^-1) to
mag/arcsec^2.
"""
function photon_radiance_to_mag(b::Bandpass, radiance)
    per_arcsec2 = uconvert(u"s^-1*m^-2", radiance * (1.0u"arcsec")^2 / u"sr")
    return _flux_to_mag(uconvert(NoUnits, per_arcsec2 / zero_point_photon_flux(b)))
end

"""
    spectral_radiance_to_mag(band, radiance)

Convert a spectral radiance (e.g. W m^-2 m^-1 sr^-1) that is flat across
the band to mag/arcsec^2.
"""
function spectral_radiance_to_mag(b::Bandpass, radiance)
    per_arcsec2 = radiance * (1.0u"arcsec")^2 / u"sr"
    return _flux_to_mag(uconvert(NoUnits, per_arcsec2 / b.zero_point))
end

# ---------------------------------------------------------------------------
# Twilight (Patat et al. 2006 V-band zenith fit; Robinson Eq. 2.84)
# ---------------------------------------------------------------------------
"""
    twilight_zenith_mag(solar_zenith_angle)

V-band zenith twilight brightness (mag/arcsec^2). Uses Patat et al. (2006)
for 95-105 deg, and the 95-deg value for smaller angles (as Robinson does).

CORRECTION: for angles above 105 deg, Robinson sets this to 22 mag/arcsec^2,
which is the whole dark-sky brightness and so counts airglow, zodiacal light
and starlight twice. Here the 105-deg flux fades linearly to zero at 108 deg
(astronomical night).
"""
function twilight_zenith_mag(solar_zenith_angle)
    g = ustrip(u"°", to_unit(u"°", solar_zenith_angle))
    patat(x) = 11.84 + 1.518 * (x - 95) - 0.057 * (x - 95)^2
    g <= 95 && return patat(95.0)
    g <= 105 && return patat(g)
    g >= 108 && return Inf
    return _flux_to_mag(_mag_to_flux(patat(105.0)) * (108 - g) / 3)
end

# ---------------------------------------------------------------------------
# Scattered moonlight (Daniels 1977, via Robinson App. B.2-B.3, Eq. 2.88-2.89)
# ---------------------------------------------------------------------------
# Lunar phase factor vs. Sun-Moon-observer angle (Robinson Table B.2).
# CORRECTION: the 120-deg entry is 0.077 in the thesis, which breaks the
# monotonic fall-off (0.041 at 110 deg, 0.017 at 130 deg) and is taken to be
# a typo. Replaced with the geometric-mean interpolation, 0.026.
const LUNAR_PHASE_ANGLE_DEG = collect(0.0:10.0:180.0)
const LUNAR_PHASE_FACTOR = [
    1.00, 0.809, 0.685, 0.483, 0.377, 0.288, 0.225, 0.172, 0.127, 0.089,
    0.061, 0.041, 0.026, 0.017, 0.009, 0.004, 0.001, 0.0, 0.0
]

# Rayleigh-scattered moonlight (Robinson Table B.3), W m^-2 m^-1 sr^-1
# (= 1e-10 W cm^-2 um^-1 sr^-1, Robinson's text units), indexed
# [z_obs, delta_az, z_moon].
const MOON_Z_OBS_DEG = collect(0.0:10.0:80.0)
const MOON_DAZ_DEG = [0.0, 45.0, 90.0, 135.0, 180.0]
const MOON_Z_MOON_DEG = [0.0, 30.0, 60.0, 75.0]
const MOON_RAYLEIGH = let
    raw = [
        [[22, 19, 13, 10], [22, 19, 13, 10], [22, 19, 13, 10], [22, 19, 13, 10], [22, 19, 13, 10]],
        [[22, 21, 15, 11], [22, 20, 14, 11], [22, 19, 13, 10], [22, 18, 12, 9.7], [22, 18, 12, 9.6]],
        [[22, 23, 18, 13], [22, 22, 16, 12], [22, 19, 14, 10], [22, 17, 12, 9.9], [22, 17, 12, 10]],
        [[22, 25, 21, 16], [22, 23, 18, 14], [22, 20, 14, 11], [22, 17, 12, 11], [22, 16, 12, 11]],
        [[23, 28, 25, 20], [23, 25, 21, 17], [23, 21, 16, 12], [23, 17, 14, 13], [23, 16, 14, 14]],
        [[24, 31, 31, 25], [24, 28, 26, 20], [24, 22, 18, 15], [24, 18, 17, 16], [24, 18, 18, 18]],
        [[27, 37, 39, 33], [27, 33, 32, 26], [27, 25, 22, 18], [27, 22, 22, 21], [27, 22, 25, 26]],
        [[34, 47, 54, 48], [34, 41, 43, 37], [34, 33, 29, 25], [34, 30, 33, 32], [34, 31, 40, 40]],
        [[55, 72, 89, 82], [55, 65, 71, 63], [55, 54, 50, 43], [55, 54, 61, 58], [58, 58, 76, 75]],
    ]
    [Float64(raw[i][j][k]) for i in 1:9, j in 1:5, k in 1:4]
end

"""
    moonlight_spectral_radiance(sky, los_zenith, los_azimuth, moon_zenith,
                                moon_azimuth, moon_separation, lunar_phase)

Scattered-moonlight spectral radiance at 550 nm (W m^-2 m^-1 sr^-1),
Robinson Eq. 2.88-2.89: `f(phase) * (F_rs + F_ms)`.

Robinson's spectral factor f_corr (Eq. 2.90) converts to a bolometric
solar-weighted signal. It is not needed here, because the radiance at 550 nm
is applied directly to the V band.

For Moon zenith angles of 75-90 deg (beyond the table) the 75-deg value fades
linearly to zero at the horizon (a modelling choice). The radiance is zero
with the Moon below the horizon.
"""
function moonlight_spectral_radiance(
    sky::SkyConditions,
    los_zenith,
    los_azimuth,
    moon_zenith,
    moon_azimuth,
    moon_separation,
    lunar_phase
)
    zm = ustrip(u"°", moon_zenith)
    zm >= 90 && return 0.0u"W/m^3/sr"
    zo = ustrip(u"°", los_zenith)
    daz = ustrip(u"°", abs(rem2pi(ustrip(u"rad", los_azimuth - moon_azimuth),
                                  RoundNearest)) * u"rad")
    horizon_fade = zm > 75 ? (90 - zm) / 15 : 1.0
    f_rs = _interp3(MOON_Z_OBS_DEG, MOON_DAZ_DEG, MOON_Z_MOON_DEG,
                    MOON_RAYLEIGH, zo, daz, min(zm, 75.0)) * horizon_fade

    psi = ustrip(u"rad", moon_separation)
    mie = sky.mie_a1 * (exp(-psi / ustrip(u"rad", sky.mie_psi1)) +
                        sky.mie_a2 * exp(-(pi - psi) / ustrip(u"rad", sky.mie_psi2)))
    phase = _interp1(LUNAR_PHASE_ANGLE_DEG, LUNAR_PHASE_FACTOR,
                     ustrip(u"°", lunar_phase))
    return phase * (1 + mie) * f_rs * u"W/m^3/sr"
end

# ---------------------------------------------------------------------------
# Zodiacal light (Roach 1972, via Robinson App. B.4)
# ---------------------------------------------------------------------------
# S10 vs. |ecliptic latitude| (rows, 0-80 deg) and elongation in ecliptic
# longitude from the Sun (columns, ascending). Robinson's table lists
# longitude descending from 180 deg; it is reversed here. "Longitude" must be
# relative to the Sun (the (0, 0) entry, 1.63e8 S10, is the Sun itself).
# Note: these values are ~1.2-1.6x brighter than Leinert et al. (1998)
# Table 16 at elongations of 90-180 deg.
const ZODI_LAT_DEG = collect(0.0:10.0:80.0)
const ZODI_DLON_DEG = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0,
                       50.0, 55.0, 60.0, 65.0, 80.0, 100.0, 120.0, 140.0, 160.0,
                       180.0]
const ZODI_S10 = let
    raw = [
        [258, 211, 207, 239, 277, 365, 535, 630, 756, 939, 1190, 1490, 2010, 2940, 4660, 7690, 15100, 36500, 176000, 163000000],
        [212, 194, 185, 217, 247, 312, 418, 455, 512, 603, 696, 825, 1150, 1550, 1820, 2140, 2760, 2720, 5630, 19900],
        [183, 174, 168, 196, 220, 258, 330, 339, 358, 403, 442, 512, 635, 800, 932, 1070, 1120, 1390, 1700, 2290],
        [159, 153, 152, 177, 196, 219, 258, 270, 282, 290, 304, 331, 363, 417, 491, 542, 592, 655, 724, 794],
        [141, 137, 137, 161, 175, 190, 204, 212, 229, 227, 233, 240, 224, 241, 246, 252, 265, 290, 315, 403],
        [127, 127, 128, 146, 156, 166, 165, 166, 183, 185, 189, 186, 171, 180, 183, 186, 190, 199, 209, 252],
        [117, 120, 120, 132, 139, 146, 137, 137, 147, 149, 150, 149, 137, 141, 144, 145, 145, 145, 146, 150],
        [110, 112, 112, 120, 123, 127, 118, 120, 124, 124, 124, 126, 118, 120, 121, 121, 121, 121, 121, 121],
        [103, 105, 105, 108, 111, 111, 106, 107, 107, 108, 107, 111, 107, 106, 108, 108, 108, 108, 108, 108],
    ]
    # Each raw row runs 180 -> 0 deg; reverse to ascending longitude.
    [Float64(reverse(raw[i])[j]) for i in 1:9, j in 1:20]
end

"""
    ecliptic_offsets(direction, sun_direction, et)

`(|ecliptic latitude|, |longitude - solar longitude|)` of `direction`, both
in degrees, the second folded into [0, 180].
"""
function ecliptic_offsets(direction, sun_direction, et)
    rot = pxform("J2000", "ECLIPJ2000", ephemeris_time(et))
    d = rot * collect(direction)
    s = rot * collect(sun_direction)
    beta = asind(clamp(d[3] / norm(d), -1.0, 1.0))
    dlon = abs(rad2deg(rem2pi(atan(d[2], d[1]) - atan(s[2], s[1]), RoundNearest)))
    return abs(beta), dlon
end

"""
    zodiacal_s10(direction, sun_direction, et)

Zodiacal-light surface brightness in S10 (outside the atmosphere).
"""
function zodiacal_s10(direction, sun_direction, et)
    beta, dlon = ecliptic_offsets(direction, sun_direction, et)
    return _interp2(ZODI_LAT_DEG, ZODI_DLON_DEG, ZODI_S10, beta, dlon)
end

# ---------------------------------------------------------------------------
# Total background
# ---------------------------------------------------------------------------
"""
    sky_background(band, sky, site, et, direction; isl_s10=0.0)

V-band sky surface brightness (mag/arcsec^2) along unit `direction` from
`site` at `et`, by component.

# Arguments
- `band::Bandpass`, `sky::SkyConditions`, `site::GroundSite`, `et`.
- `direction`: Unit line-of-sight vector (J2000).
- `isl_s10`: Integrated starlight above the atmosphere, in S10, from stars
  fainter than those drawn explicitly (see `integrated_starlight_s10`).

# Returns
A NamedTuple of brightnesses (mag/arcsec^2; `Inf` means none) for
`airglow`, `light_pollution`, `twilight`, `moonlight`, `zodiacal`,
`starlight` and `total`, plus `zenith_angle`, `airmass` and `transmission`.
"""
function sky_background(
    b::Bandpass,
    sky::SkyConditions,
    site::GroundSite,
    et,
    direction::AbstractVector;
    isl_s10::Real = 0.0
)
    observer = site_state(site, et)
    sun = spice_state("SUN", et)
    moon = spice_state("MOON", et)

    los = horizon_coordinates(site, et, direction)
    airmass = relative_airmass(los.zenith)
    transmission = atmospheric_transmission(sky, los.zenith)

    # Atmospheric emission/scattering: zenith value x relative airmass.
    airglow = photon_radiance_to_mag(b, airglow_zenith_photon_radiance(sky, b)) -
              2.5 * log10(airmass)
    pollution = sky.light_pollution === nothing ? Inf :
                sky.light_pollution - 2.5 * log10(airmass)
    sun_zenith = horizon_coordinates(site, et, line_of_sight(sun, observer)).zenith
    twilight = twilight_zenith_mag(sun_zenith) - 2.5 * log10(airmass)

    moon_dir = line_of_sight(moon, observer)
    moon_hz = horizon_coordinates(site, et, moon_dir)
    lunar_phase = solar_phase_angle(moon, sun, observer) * u"rad"
    separation = acos(clamp(dot(normalize(SVector{3}(direction)), moon_dir),
                            -1.0, 1.0)) * u"rad"
    moonlight = spectral_radiance_to_mag(
        b,
        moonlight_spectral_radiance(sky, los.zenith, los.azimuth, moon_hz.zenith,
                                    moon_hz.azimuth, separation, lunar_phase)
    )

    # Outside the atmosphere: attenuated by extinction.
    sun_dir = line_of_sight(sun, observer)
    zodiacal = s10_to_mag(zodiacal_s10(direction, sun_dir, et) * transmission)
    starlight = s10_to_mag(isl_s10 * transmission)

    components = (
        airglow = airglow,
        light_pollution = pollution,
        twilight = twilight,
        moonlight = moonlight,
        zodiacal = zodiacal,
        starlight = starlight
    )
    total = _flux_to_mag(sum(_mag_to_flux, values(components)))
    return merge(components, (
        total = total,
        zenith_angle = los.zenith,
        airmass = airmass,
        transmission = transmission,
        moon_separation = separation,
        moon_zenith = moon_hz.zenith,
        lunar_phase = lunar_phase
    ))
end
