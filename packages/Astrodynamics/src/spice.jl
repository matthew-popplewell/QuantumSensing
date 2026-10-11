"""
    furnsh_kernels(kernels::Vector{String} = [])

Download and furnsh SPICE kernels.

# Arguments
- `kernels::Vector{String}`: A vector of URLs or file paths to the SPICE kernels to be downloaded and furnshed.
"""
function furnsh_kernels(kernels::Vector{String} = [
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"
])
    # Download and furnsh SPICE kernels
    kernel_dir = joinpath(@__DIR__, "../../../data/spice")
    if !isdir(kernel_dir)
        mkpath(kernel_dir)
    end

    for kernel in kernels
        kernel_file = joinpath(kernel_dir, basename(kernel))
        if !isfile(kernel_file)
            println("Downloading $kernel...")
            download(kernel, kernel_file)
        end
        SPICE.furnsh(kernel_file)
    end
end

"""
    build_spline(tspan::Tuple{Float64, Float64}, body_id::String)

Build a spline for the specified celestial body over the given time span.

# Arguments
- `tspan::Tuple{Float64, Float64}`: A tuple specifying the start and end times for the spline in ephemeris time (ET).
- `body_id::String`: The NAIF ID of the celestial body for which to build the spline.
- `step_size::Float64`: The time step size for the spline (default is 60.0 seconds).
- `reference_frame::String`: The reference frame for the state vectors (default is "J2000").
- `observer::String`: The observer body for the state vectors (default is "EARTH").
- `light_time_correction::String`: The light time correction to apply (default is "NONE").
"""
function build_spline(
    tspan::Tuple{Float64, Float64},
    body_id::String,
    step_size::Float64 = 60.0,
    reference_frame::String = "J2000",
    observer::String = "EARTH",
    light_time_correction::String = "NONE"
)
    tspan = collect(tspan[1]-0.1:step_size:tspan[2]+0.1)
    states = first.(SPICE.spkezr(body_id, tspan, reference_frame, light_time_correction, observer))
    return CubicSpline(tspan, states)
end
