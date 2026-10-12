"""
    solar_phase_angle(target, sun, observer)
Compute the solar phase angle for a target object given the positions of the sun and an observer.

# Arguments
- `target::SVector{6}`: The state vector of the target object [m, m/s].
- `sun::SVector{6}`: The state vector of the sun [m, m/s].
- `observer::SVector{6}`: The state vector of the observer [m, m/s].
"""
function relative_directions(
    target::SVector{6},
    sun::SVector{6},
    observer::SVector{6}
)
    r_target = target[1:3]
    r_sun = sun[1:3]
    r_observer = observer[1:3]
    object2sun = (r_sun - r_target) / norm(r_sun - r_target)
    object2observer = (r_observer - r_target) / norm(r_observer - r_target)
    return SVector{3}(object2sun), SVector{3}(object2observer)
end

"""
    solar_phase_angle(object2sun, object2observer)
Compute the solar phase angle given the unit vectors from the target to the sun and from the target
to the observer.

# Arguments
- `object2sun::SVector{3}`: The unit vector from the target to the sun.
- `object2observer::SVector{3}`: The unit vector from the target to the observer.
"""
function solar_phase_angle(object2sun::SVector{3}, object2observer::SVector{3})
    return acos(clamp(dot(object2sun, object2observer), -1.0, 1.0))
end

"""
    solar_apparent_radius(target, sun)
Compute the apparent radius of the sun as seen from the target object.

# Arguments
- `target::SVector{6}`: The state vector of the target object [m, m/s].
- `sun::SVector{6}`: The state vector of the sun [m, m/s].
"""
function solar_apparent_radius(target::SVector{6}, sun::SVector{6})
    _, sun_radii_km = bodvrd("SUN", "RADII", 3)
    sun_radius = sun_radii_km[1] * 1e3  # SPICE returns km; convert to m
    distance = norm(sun[1:3] - target[1:3])
    return asin(sun_radius / distance)
end

"""
    lambertian_extended_source(target, sun, observer, radius, coeff_diffusion)
Compute the exitance of a Lambertian extended source

# Arguments
- `target::SVector{6}`: The state vector of the target object [m, m/s].
- `sun::SVector{6}`: The state vector of the sun [m, m/s].
- `observer::SVector{6}`: The state vector of the observer [m, m/s].
- `radius::Float64`: The radius of the target object [m].
- `coeff_diffusion::Float64`: The diffusion coefficient of the target object's surface.
"""
function lambertian_extended_source(
    target::SVector{6},
    sun::SVector{6},
    observer::SVector{6},
    radius::Float64,
    coeff_diffusion::Float64
)
    # Relative directions
    object2sun, object2observer = relative_directions(target, sun, observer)

    # Angles
    phase_angle = solar_phase_angle(object2sun, object2observer)
    apparent_radius = solar_apparent_radius(target, sun)

    exitance = 1/6 * coeff_diffusion/pi * radius^2 * cos(apparent_radius) *
        (5-cos(apparent_radius)) * (sin(phase_angle)+ (pi-phase_angle)*cos(phase_angle))

    return exitance
end