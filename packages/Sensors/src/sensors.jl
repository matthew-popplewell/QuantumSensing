"""
    CCD(aperature_area, effective_aperture_area, focal_length, pixel_scale, quantum_efficiency, dark_current, read_noise, exposure_time, gain)
CCD sensor parameters.

# Fields
- `aperature_area::Float64`: The area of the CCD's aperature.
- `effective_aperture_area::Float64`: The effective area of the CCD's aperature.
- `focal_length::Float64`: The focal length of the CCD.
- `pixel_scale::Float64`: The pixel scale of the CCD.
- `quantum_efficiency::Float64`: The quantum efficiency of the CCD (between 0 and 1).
- `dark_current::Float64`: The dark current of the CCD in electrons per second per pixel.
- `read_noise::Float64`: The read noise of the CCD in electrons per pixel.
- `exposure_time::Float64`: The exposure time of the CCD.
- `gain::Float64`: The gain of the CCD in electrons per ADU (analog-to-digital unit).
- `effective_wavelength::Float64`: The effective wavelength of the CCD.
- `fwhm::Float64`: The Full Width at Half Maximum of the PSF in arcseconds.
- `boresight_ra::Float64`: The boresight right ascension of the CCD in radians.
- `boresight_dec::Float64`: The boresight declination of the CCD in radians.
"""
struct CCD
    aperture_area::Quantity{Float64, Unitful.𝐋^2}
    effective_aperture_area::Quantity{Float64, Unitful.𝐋^2}
    focal_length::Quantity{Float64, Unitful.𝐋}
    pixel_scale::Quantity{Float64, Unitful.NoDims}
    pixel_pitch::Quantity{Float64, Unitful.𝐋}
    quantum_efficiency::Float64
    dark_current::Quantity{Float64, Unitful.𝐓^-1}
    read_noise::Float64
    exposure_time::Quantity{Float64, Unitful.𝐓}
    gain::Float64
    effective_wavelength::Quantity{Float64, Unitful.𝐋}
    FWHM::Quantity{Float64, Unitful.NoDims}
    boresight_ra::Quantity{Float64, Unitful.NoDims}
    boresight_dec::Quantity{Float64, Unitful.NoDims}
end

"""
    compute_start_pixel(target, observer, boresight_ra_dec, focal_length, pixel_pitch, image_size)

Projects a target's 3D state vector onto the CCD pixel grid given telescope pointing.
"""
function compute_start_pixel(
    ccd::CCD,
    target::AbstractVector{<:Quantity},
    observer::AbstractVector{<:Quantity},
    image_size::Tuple{Int, Int}
)
    # Strip
    alpha_0 = ccd.boresight_ra
    delta_0 = ccd.boresight_dec
    focal_length = ustrip(u"m", ccd.focal_length)
    pixel_pitch = ustrip(u"m", ccd.pixel_pitch)

    # Convert positions to meters for uniform calculation
    target_m = ustrip.(u"m", target[1:3])
    observer_m = ustrip.(u"m", observer[1:3])
    
    rho = target_m - observer_m
    rho_norm = norm(rho)

    # Target unit vector (dimensionless components)
    ux = ustrip(rho[1] / rho_norm)
    uy = ustrip(rho[2] / rho_norm)
    uz = ustrip(rho[3] / rho_norm)

    delta = asin(uz)
    alpha = atan(uy, ux)

    # Gnomonic (Tangent Plane) Projection
    cos_d, sin_d = cos(delta), sin(delta)
    cos_d0, sin_d0 = cos(delta_0), sin(delta_0)
    cos_d_ra, sin_d_ra = cos(alpha - alpha_0), sin(alpha - alpha_0)

    denominator = sin_d0 * sin_d + cos_d0 * cos_d * cos_d_ra
    xi = (cos_d * sin_d_ra) / denominator
    eta = (cos_d0 * sin_d - sin_d0 * cos_d * cos_d_ra) / denominator

    # Convert to focal plane meters
    x_focal = focal_length * xi
    y_focal = focal_length * eta

    # Convert to pixel coordinates relative to chip center
    nrows, ncols = image_size
    cx = ncols / 2.0
    cy = nrows / 2.0

    x_pixel = cx + ustrip(x_focal / pixel_pitch)
    y_pixel = cy + ustrip(y_focal / pixel_pitch)

    return (x_pixel, y_pixel)
end

"""
    compute_pixel_velocity(ccd, target, observer, tracking_mode)

Computes the apparent pixel velocity (pixels/second) of a target object accounting 
for telescope mount tracking modes (:sidereal or :geo_track).

# Arguments
- `ccd::CCD`: The CCD sensor parameters.
- `target::SVector{6, Float64}`: The state vector (position and velocity) of the target object.
- `observer::SVector{6, Float64}`: The state vector (position and velocity) of the observer.
- `tracking_mode::Symbol`: The tracking mode of the telescope mount. Can be `:sidereal` or `:geo_track`.
"""
function compute_pixel_velocity(
    ccd::CCD,
    target::AbstractVector{<:Quantity},       
    observer::AbstractVector{<:Quantity},
    tracking_mode::Symbol = :sidereal  
)
    # Standardize position to meters and velocity to m/s
    rho = ustrip.(u"m", target[1:3] - observer[1:3])
    v_rel = ustrip.(u"m/s", target[4:6] - observer[4:6])
    range_dist = norm(rho)

    # 2. Target inertial apparent angular velocity vector [rad/s]
    omega_target = cross(rho, v_rel) / (range_dist^2)

    # 3. Telescope tracking rate
    omega_earth = 7.292115e-5u"rad/s"
    if tracking_mode == :sidereal
        omega_scope = SVector(0.0u"rad/s", 0.0u"rad/s", omega_earth)
    elseif tracking_mode == :geo_track
        omega_scope = omega_target
    else
        error("Unknown tracking mode: $tracking_mode. Use :sidereal or :geo_track.")
    end

    # 4. Effective angular velocity [rad/s]
    omega_eff = omega_target - omega_scope

    # 5. Convert to focal plane linear velocity and then pixel velocity
    v_focal_x = omega_eff[2] * ccd.focal_length
    v_focal_y = -omega_eff[1] * ccd.focal_length

    vx_pix = ustrip(v_focal_x / ccd.pixel_pitch)
    vy_pix = ustrip(v_focal_y / ccd.pixel_pitch)

    return (vx_pix, vy_pix)
end

"""
    simulate_ccd_signal(ccd, target, sun, observer, radius, coeff_diffusion)
Simulate the signal received by a CCD sensor from a target object.

# Arguments
- `ccd::CCD`: The CCD sensor parameters.
- `target::SVector{6, Float64}`: The state vector (position and velocity) of the target object.
- `sun::SVector{6, Float64}`: The state vector (position and velocity) of the sun.
- `observer::SVector{6, Float64}`: The state vector (position and velocity) of the observer.
- `radius::Float64`: The radius of the target object.
- `coeff_diffusion::Float64`: The diffusion coefficient of the target object's surface.
"""
function simulate_ccd_signal(
    ccd::CCD,
    target::SVector{6},
    sun::SVector{6},
    observer::SVector{6},
    radius::Float64,
    coeff_diffusion::Float64
)
    # Extract SI base values safely using ustrip
    effective_aperture_area = ustrip(u"m^2", ccd.effective_aperture_area)
    exposure_time = ustrip(u"s", ccd.exposure_time)
    effective_wavelength = ustrip(u"m", ccd.effective_wavelength)

    r_sun = ustrip.(u"m", sun[1:3])
    r_target = ustrip.(u"m", target[1:3])
    r_observer = ustrip.(u"m", observer[1:3])
    
    au_m = ustrip(u"m", AU)
    solar_flux = ustrip(u"W/m^2", SOLAR_CONSTANT)
    planck = ustrip(u"J*s", PLANCKS_CONSTANT)
    speed_light = ustrip(u"m/s", SPEED_LIGHT)

    # Reflectance of the target object using Lambertian model
    reflectance = lambertian_extended_source(target, sun, observer, radius, coeff_diffusion)
    
    # Irradiance at the CCD
    solar_irradiance = solar_flux * (au_m / norm(r_sun - r_target))^2
    target_irradiance = reflectance * solar_irradiance / norm(r_observer - r_target)^2

    # Photon flux at the CCD
    photon_energy = planck * speed_light / effective_wavelength
    photon_flux = target_irradiance * effective_aperture_area / photon_energy # whole aperature
    
    # Deficiencies and exposure time
    signal = photon_flux * ccd.quantum_efficiency * exposure_time

    return signal
end

function simulate_ccd_image(
    ccd::CCD,
    target::SVector{6},
    sun::SVector{6},
    observer::SVector{6},
    radius::Float64,
    coeff_diffusion::Float64,
    image_size::Tuple{Int, Int} = (1024, 1024)
)
    nrows, ncols = image_size
    image = zeros(Float64, nrows, ncols)

    exposure_time = ustrip(u"s", ccd.exposure_time)
    dark_current = ustrip(u"s^-1", ccd.dark_current)
    boresight_ra = ustrip(u"rad", ccd.boresight_ra)
    boresight_dec = ustrip(u"rad", ccd.boresight_dec)
    fwhm = ustrip(ccd.FWHM)
    pixel_scale = ustrip(ccd.pixel_scale)

    # Signal
    target_signal = simulate_ccd_signal(
        ccd, target, sun, observer, radius, coeff_diffusion
    )

    # Background signal and dark current
    background_signal = 0.0 # Placeholder for background signal calculation
    noise_lambda = exposure_time * (background_signal + dark_current)
    noise = rand(Poisson(noise_lambda), (nrows, ncols))
    image .+= noise
    
    # Sub-stepped spatial projection
    pixel_velocity = compute_pixel_velocity(ccd, target, observer)
    start_pixel = compute_start_pixel(
        ccd, target, observer, image_size
    )
    fwhm_pixels = fwhm / pixel_scale
    sigma = fwhm_pixels / (2.0 * sqrt(2.0 * log(2.0)))
    displacement_magnitude = sqrt(pixel_velocity[1]^2 + pixel_velocity[2]^2) * exposure_time
    num_substeps = max(10, Int(ceil(displacement_magnitude * 4.0)))
    dt = exposure_time / num_substeps
    signal_per_substep = target_signal / num_substeps
    box_radius = ceil(Int, 3.0 * sigma)

    for step in 0:(num_substeps-1)
        t = (step + 0.5) * dt
        xc = start_pixel[1] + pixel_velocity[1] * t
        yc = start_pixel[2] + pixel_velocity[2] * t

        min_r = max(1, floor(Int, yc - box_radius))
        max_r = min(nrows, ceil(Int, yc + box_radius))
        min_c = max(1, floor(Int, xc - box_radius))
        max_c = min(ncols, ceil(Int, xc + box_radius))

        local_sum = 0.0
        local_weights = zeros(Float64, max_r - min_r + 1, max_c - min_c + 1)

        for r in min_r:max_r, c in min_c:max_c
            dx = c - xc
            dy = r - yc
            w = exp(-(dx^2 + dy^2) / (2.0 * sigma^2))
            local_weights[r - min_r + 1, c - min_c + 1] = w
            local_sum += w
        end

        if local_sum > 0
            for r in min_r:max_r, c in min_c:max_c
                image[r, c] += signal_per_substep * (local_weights[r - min_r + 1, c - min_c + 1] / local_sum)
            end
        end
    end

    # Inject Poisson shot noise (on true electrons) and Gaussian read noise per pixel
    noisy_image = zeros(Float64, nrows, ncols)
    for r in 1:nrows, c in 1:ncols
        lambda = max(0.0, image[r, c])
        actual_electrons = rand(Poisson(lambda))
        read_noise = rand(Normal(0.0, ccd.read_noise))
        noisy_image[r, c] = actual_electrons + read_noise
    end

    # Electron counts to final digital numbers
    adu_image = round.(Int, noisy_image ./ ccd.gain)

    return adu_image

end