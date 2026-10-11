using Sensors
using Unitful
using LinearAlgebra
using StaticArrays
using SPICE
using Astrodynamics

# Furnsh
furnsh_kernels()

D_tel = 0.5u"m"
ap_area = (pi / 4.0) * D_tel^2

my_ccd = CCD(
    ap_area, 
    ap_area * 0.85, 
    4.0u"m", 
    0.5, 
    9.0u"μm", 
    0.80, 
    0.02u"s^-1",
    3.5, 
    2.0u"s", 
    2.3, 
    550.0u"nm", 
    1.5, 
    0.0u"rad", 
    0.0u"rad"
)

println("Defining state vectors in kilometers and km/s...")

# Observer on Earth surface defined in kilometers and km/s
observer_state = SVector(6378.137u"km", 0.0u"km", 0.0u"km", 0.0u"km/s", 0.4651u"km/s", 0.0u"km/s")

# GEO satellite defined in kilometers and km/s
target_state = SVector(42164.0u"km", 1.0u"km", 0.5u"km", -0.01u"km/s", 3.07u"km/s", 0.1u"km/s")

# Sun defined in AU or km
sun_state = SVector(AU, 0.0u"km", 0.0u"km", 0.0u"km/s", 29.78u"km/s", 0.0u"km/s")

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