module Sensors

using LinearAlgebra
using Unitful
using StaticArrays
using Random
using Distributions
using SPICE

# Unitful does not ship an arcsecond unit; define it here.
@unit arcsec "arcsec" Arcsecond (1 / 3600)u"°" false
Unitful.register(Sensors)

function __init__()
    Unitful.register(Sensors)
end

# Constants (Unitful's own `Unitful.h` and `Unitful.c0` are used for the
# Planck constant and speed of light).
const AU = 149_597_870_700.0u"m"                 # IAU 2012 Resolution B2
const EARTH_GM = 398_600.4418u"km^3/s^2"         # WGS 84
const EARTH_ROTATION_RATE = 7.292115e-5u"rad/s"  # WGS 84
const GEO_RADIUS = uconvert(
    u"km", cbrt(EARTH_GM / uconvert(u"s^-1", EARTH_ROTATION_RATE)^2)
)

include("units.jl")
include("geometry.jl")
include("constraints.jl")
include("brdf.jl")
include("sensors.jl")
include("background.jl")
include("stars.jl")
include("signals.jl")

# Bessell V. Pivot wavelength, effective width and zero point (3630.2 Jy, Vega)
# from the SVO Filter Profile Service, Generic/Bessell.V:
# http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Generic/Bessell.V
# Solar V magnitude (Vega) from Willmer (2018), ApJS 236, 47, Table 3.
const BESSELL_V = Bandpass(
    wavelength = 5477.6973781164u"Å",
    width = 893.09998496416u"Å",
    zero_point = 3630.2172842325e-26u"W/m^2/Hz",
    sun_magnitude = -26.76
)

# Constants
export AU
export EARTH_GM
export EARTH_ROTATION_RATE
export GEO_RADIUS
export BESSELL_V

# units.jl
export CartesianState
export state_vector
export drift

# geometry.jl
export ephemeris_time
export spice_state
export GroundSite
export site_state
export site_zenith
export geostationary_state
export line_of_sight
export radec

# constraints.jl
export elevation
export in_earth_shadow
export illumination_fraction

# brdf.jl
export relative_directions
export solar_phase_angle
export solar_apparent_radius
export lambertian_extended_source

# sensors.jl
export Telescope
export Detector
export Bandpass
export Sensor
export Observation
export collecting_area
export pixel_scale
export field_of_view
export point_at
export boresight_at
export target_pixel
export sky_pixel
export pixel_velocity

# signals.jl
export photon_energy
export zero_point_photon_flux
export apparent_magnitude
export electron_rate
export sky_electron_rate
export simulate_image

# background.jl
export SkyConditions
export relative_airmass
export atmospheric_transmission
export horizon_coordinates
export twilight_zenith_mag
export zodiacal_s10
export s10_to_mag
export sky_background

# stars.jl
export SyntheticStar
export gaia_g_minus_v
export galactic_coordinates
export star_counts
export integrated_starlight_s10
export synthetic_star_field

end # module Sensors
