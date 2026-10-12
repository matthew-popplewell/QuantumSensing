using Sensors
using LinearAlgebra
using Statistics
using StaticArrays
using SPICE
using Astrodynamics

# Furnsh
furnsh_kernels()

# All quantities are in SI units: m, s, rad.
ARCSEC = pi / (180.0 * 3600.0)  # [rad/arcsec]

D_tel = 0.5  # [m]
ap_area = (pi / 4.0) * D_tel^2  # [m^2]

my_ccd = CCD(
    ap_area,           # Aperture area [m^2]
    ap_area * 0.85,    # Effective aperture area [m^2]
    4.0,               # Focal length [m]
    0.5 * ARCSEC,      # Pixel scale [rad]
    9.0e-6,            # Pixel pitch [m]
    0.80,              # Quantum efficiency
    0.02,              # Dark current [electrons/s/pixel]
    3.5,               # Read noise [electrons/pixel]
    2.0,               # Exposure time [s]
    2.3,               # Gain [electrons/ADU]
    550.0e-9,          # Effective wavelength [m]
    1.5 * ARCSEC,      # PSF FWHM [rad]
    0.0,               # Boresight RA [rad]
    0.0                # Boresight Dec [rad]
)

println("Defining state vectors in meters and m/s...")

# [x, y, z, vx, vy, vz] in [m, m, m, m/s, m/s, m/s]
observer_state = SVector{6, Float64}(
    6378.137e3, 0.0, 0.0,
    0.0, 465.1, 0.0
)

target_state = SVector{6, Float64}(
    42164.0e3, 1.0e3, 0.5e3,
    -10.0, 3070.0, 100.0
)

sun_state = SVector{6, Float64}(
    149_597_870.7e3, 0.0, 0.0,
    0.0, 29.78e3, 0.0
)

target_radius = 2.5       # [m]
coeff_diffusion = 0.4     # Albedo

println("Running CCD image simulation...")
simulated_frame = simulate_ccd_image(
    my_ccd,
    target_state,
    sun_state,
    observer_state,
    target_radius,
    coeff_diffusion,
    (512, 512),
)

# Print diagnostic statistics about the resulting simulated image frame
println("Simulation complete!")
println("Image Dimensions: $(size(simulated_frame))")
println("Min Pixel Value (ADU): $(minimum(simulated_frame))")
println("Max Pixel Value (ADU): $(maximum(simulated_frame))")
println("Mean Pixel Value (ADU): $(round(mean(simulated_frame), digits=2))")

# Find the peak center pixel where the target landed
peak_val, peak_idx = findmax(simulated_frame)
println("Target Peak Detected at Pixel Coordinate: $peak_idx with intensity $peak_val ADU")