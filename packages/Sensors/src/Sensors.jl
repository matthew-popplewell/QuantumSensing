module Sensors

using LinearAlgebra
using StaticArrays
using Random
using Distributions
using SPICE

# Physical constants (SI units)
PLANCKS_CONSTANT = 6.62607015e-34  # [J s]
SPEED_LIGHT = 299792458.0          # [m/s]
MAG_SUN = -26.74                   # Apparent magnitude of the Sun
AU = 149597870700.0                # [m]
SOLAR_CONSTANT = 1361.0            # [W/m^2]

include("brdf.jl")
include("signals.jl")
include("sensors.jl")

# Constants
export PLANCKS_CONSTANT
export SPEED_LIGHT
export MAG_SUN
export AU
export SOLAR_CONSTANT

# brdf.jl
export relative_directions
export solar_phase_angle
export solar_apparent_radius
export lambertian_extended_source

# signals.jl

# sensors.jl
export CCD
export compute_start_pixel
export compute_pixel_velocity
export simulate_ccd_signal
export simulate_ccd_image

end # module Sensors
