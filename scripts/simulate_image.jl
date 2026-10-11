#=
Simulate a ground-based image of a GEO satellite.

Station:   Purdue Optical Ground Station (POGS), New Mexico, with the telescope
           and sensor parameters of Robinson (2023), "Light curve simulation
           and shape inversion for human-made space objects", MS thesis,
           Purdue University, Appendix C, Table C.1. Bessell V band.
Target:    Ideal geostationary Lambertian sphere near the site's meridian.
Epoch:     The night of 2026-10-10/11; the script picks the lowest-phase-angle
           time that is astronomically dark, above 30 deg elevation and fully
           sunlit (mid-October is GEO eclipse season, so local midnight is in
           shadow).

Sky background: Robinson (2023) Sec. 2.5.5 (airglow, light pollution,
twilight, moonlight, zodiacal light, integrated starlight), with the
corrections documented in packages/Sensors/src/background.jl. Stars:
synthetic field with Gaia DR3-calibrated counts.

Values marked [ASSUMED] are modelling choices and [INFERRED] are not in
Robinson's table; neither is a sourced POGS figure.
=#

using Sensors
using LinearAlgebra
using Unitful
using SPICE
using Statistics
using Astrodynamics

furnsh_kernels()

# ---------------------------------------------------------------------------
# Instrument: POGS, Robinson (2023) Table C.1
# ---------------------------------------------------------------------------
# Aperture 0.35560 m, f/7.2, secondary 0.1724660 m. Robinson's signal model
# has no optical losses (only the effective aperture area, QE and the
# atmosphere), so throughput = 1 to match it.
aperture = 0.35560u"m"
telescope = Telescope(
    aperture_diameter = aperture,
    focal_length = 7.2 * aperture,
    obstruction_ratio = 0.1724660 / 0.35560,
    throughput = 1.0
)

# 4096 x 4096 px of 9.009 um, dark current 3 e-/px/s, read-noise variance
# 9 e-^2, gain 1 e-/ADU (Table C.1).
# Robinson uses Krag's (2003) QE curve without tabulating it. A 4096^2 array
# of 9 um pixels matches the KAF-16803 sensor, so its datasheet QE (60% at
# 550 nm) and full well (100k e-) are used [INFERRED]:
#   https://web.archive.org/web/20190621014904/https://www.onsemi.com/pub/Collateral/KAF-16803-D.PDF
detector = Detector(
    pixel_pitch = 9.009u"μm",
    npixels = (4096, 4096),
    quantum_efficiency = 60u"percent",  # [INFERRED]
    read_noise = sqrt(9.0),
    dark_current = 3.0u"s^-1",
    gain = 1.0,
    full_well = 100_000,                # [INFERRED]
    bias = 1000.0,                      # [ASSUMED]
    bit_depth = 16                      # [ASSUMED]
)

sensor = Sensor(telescope, detector, BESSELL_V)
frame_size = detector.npixels  # full 0.82 deg x 0.82 deg frame

# ---------------------------------------------------------------------------
# Site and target
# ---------------------------------------------------------------------------
# POGS location from Robinson (2023), e.g. Fig. 2.17 caption: 32.900 N,
# -105.533 W. The thesis does not state the altitude in text; 2248.7 m above
# mean sea level is from the USGS 3DEP Elevation Point Query Service at that
# point (https://epqs.nationalmap.gov/v1/json?x=-105.533&y=32.900&units=Meters).
# The geoid-ellipsoid offset (tens of metres) is ignored.
site = GroundSite(latitude = 32.900u"°", longitude = -105.533u"°",
                  altitude = 2248.7u"m")
geo_longitude = -105.0u"°"  # [ASSUMED] synthetic GEO slot near the meridian
target_radius = 1.5u"m"     # [ASSUMED]
target_albedo = 0.2         # [ASSUMED]

# ---------------------------------------------------------------------------
# Pick the observation time
# ---------------------------------------------------------------------------
function observable_geometry(et)
    observer = site_state(site, et)
    target = geostationary_state(longitude = geo_longitude, et = et)
    sun = spice_state("SUN", et)
    return (
        et = et,
        observer = observer,
        target = target,
        sun = sun,
        sun_elevation = elevation(site, et, sun),
        target_elevation = elevation(site, et, target),
        illumination = illumination_fraction(target, sun),
        phase = solar_phase_angle(target, sun, observer) * u"rad"
    )
end

et_start = ephemeris_time("2026-10-10 19:00:00 UTC")
candidates = [observable_geometry(et_start + 300.0 * k) for k in 0:144]
usable = filter(
    g -> g.sun_elevation < -18u"°" && g.target_elevation > 30u"°" &&
         g.illumination == 1.0,
    candidates
)
isempty(usable) && error("No usable observation time found.")
geom = argmin(g -> g.phase, usable)

println("Epoch:              ", et2utc(geom.et, "ISOC", 0), " UTC")
println("Sun elevation:      ", round(u"°", geom.sun_elevation, digits = 1))
println("Target elevation:   ", round(u"°", geom.target_elevation, digits = 1))
println("Phase angle:        ", round(u"°", geom.phase, digits = 1))
println("Range:              ",
        round(u"km", norm(geom.target.r - geom.observer.r), digits = 0))

# ---------------------------------------------------------------------------
# Sky conditions
# ---------------------------------------------------------------------------
sky = SkyConditions(
    extinction = 0.11,        # Robinson Table B.1 at 550 nm (mag/airmass)
    # Artificial zenith brightness, V mag/arcsec^2. Robinson interpolates the
    # Falchi et al. (2016) World Atlas; look up POGS there (e.g.
    # lightpollutionmap.info) and set it here. `nothing` = none. [ASSUMED]
    light_pollution = nothing,
    sodium_rayleigh = 65.0    # Leinert (1998) Table 13: 30 (summer)-100 (winter) R
)

# ---------------------------------------------------------------------------
# Simulate
# ---------------------------------------------------------------------------
println("\nPlate scale:        ", round(u"arcsec", pixel_scale(sensor), digits = 3))
println("Subframe FOV:       ",
        round.(u"°", field_of_view(sensor, frame_size), digits = 3))

"""
    aperture_snr(img, sigma_px, background_e)

SNR from Robinson Eq. 2.70 for an aperture of radius one FWHM around the
target (a circle for a point source, a capsule along a streak).
"""
function aperture_snr(img, sigma_px, background_e, read_noise, gain)
    r = 2 * sqrt(2 * log(2)) * sigma_px
    len = hypot((img.end_pixel .- img.start_pixel)...)
    npix = pi * r^2 + 2 * r * len
    signal = img.target_electrons * (1 - exp(-r^2 / (2 * sigma_px^2)))
    noise_var = npix * (background_e + read_noise^2 + gain^2 / 12)
    return signal / sqrt(signal + noise_var)
end

images = Dict{Symbol, Any}()
for tracking in (:sidereal, :rate_track)
    obs = Observation(
        boresight = point_at(geom.target, geom.observer),
        exposure_time = 10u"s",       # Robinson Table C.1
        seeing_fwhm = 1.5u"arcsec",   # Robinson Table C.1 FWHM
        tracking = tracking
    )

    img = simulate_image(
        sensor, obs, site, geom.et, geom.target;
        radius = target_radius,
        albedo = target_albedo,
        sky = sky,
        frame_size = frame_size
    )
    images[tracking] = img

    bg = img.background
    bg_e = img.background_electrons
    sigma_px = uconvert(NoUnits, obs.seeing_fwhm / pixel_scale(sensor)) /
               (2 * sqrt(2 * log(2)))
    vx, vy = pixel_velocity(sensor, obs, geom.target, geom.observer)
    rate = hypot(vx, vy) * pixel_scale(sensor) / 1u"s"

    println("\n== Tracking: $tracking ==")
    println("V magnitude (exo-atm.): ", round(img.magnitude, digits = 2),
            "   illumination: ", round(img.illumination, digits = 3))
    println("Zenith angle / airmass: ", round(u"°", bg.zenith_angle, digits = 1),
            " / ", round(bg.airmass, digits = 3),
            "   transmission: ", round(bg.transmission, digits = 3))
    println("Moon: zenith ", round(u"°", bg.moon_zenith, digits = 1),
            ", phase angle ", round(u"°", bg.lunar_phase, digits = 1),
            ", separation ", round(u"°", bg.moon_separation, digits = 1))
    println("Target electrons:       ", round(Int, img.target_electrons))
    println("Apparent rate:          ", round(u"arcsec/s", rate, digits = 2))
    println("Streak:                 ", round.(img.start_pixel, digits = 1),
            " -> ", round.(img.end_pixel, digits = 1))
    println("Stars drawn:            ", length(img.stars),
            " (V < ", round(img.star_limit, digits = 2), ")")
    println("Background budget:")
    for k in (:airglow, :light_pollution, :twilight, :moonlight, :zodiacal,
              :starlight)
        println("  ", rpad(string(k), 16),
                rpad(isfinite(bg[k]) ? string(round(bg[k], digits = 2)) * " mag/as2" : "none", 18),
                round(bg_e[k], digits = 2), " e-/px")
    end
    sky_e = sum(bg_e[k] for k in keys(bg_e) if k != :dark)
    println("  ", rpad("total sky", 16),
            rpad(string(round(bg.total, digits = 2)) * " mag/as2", 18),
            round(sky_e, digits = 2), " e-/px")
    println("  ", rpad("dark", 34), round(bg_e.dark, digits = 3), " e-/px")
    println("Aperture SNR (Eq. 2.70): ",
            round(aperture_snr(img, sigma_px, sky_e + bg_e.dark,
                               detector.read_noise, detector.gain), digits = 1))
    println("ADU min/median/max:     ", minimum(img.adu), " / ",
            median(img.adu), " / ", maximum(img.adu))
end
