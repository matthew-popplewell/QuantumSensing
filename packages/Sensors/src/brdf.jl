"""
    relative_directions(target, sun, observer)

Unit vectors from the target to the Sun and from the target to the observer.
States may be `CartesianState`s or 6-vectors of quantities.
"""
function relative_directions(target, sun, observer)
    r_target = CartesianState(target).r
    to_sun = CartesianState(sun).r - r_target
    to_observer = CartesianState(observer).r - r_target
    return (
        SVector{3, Float64}(to_sun / norm(to_sun)),
        SVector{3, Float64}(to_observer / norm(to_observer))
    )
end

"""
    solar_phase_angle(object2sun, object2observer)
    solar_phase_angle(target, sun, observer)

Sun-target-observer angle (rad, plain Float64).
"""
function solar_phase_angle(object2sun::SVector{3}, object2observer::SVector{3})
    return acos(clamp(dot(object2sun, object2observer), -1.0, 1.0))
end

function solar_phase_angle(target, sun, observer)
    return solar_phase_angle(relative_directions(target, sun, observer)...)
end

"""
    solar_apparent_radius(target, sun)

Angular radius of the Sun seen from the target (rad, plain Float64).
"""
function solar_apparent_radius(target, sun)
    sun_radius = bodvrd("SUN", "RADII")[1] * u"km"
    distance = norm(CartesianState(sun).r - CartesianState(target).r)
    return asin(uconvert(NoUnits, sun_radius / distance))
end

"""
    lambertian_extended_source(target, sun, observer, radius, albedo)

Reflection factor of a Lambertian sphere illuminated by the (extended) Sun:
the reflected irradiance at distance `d` is
`solar_irradiance * lambertian_extended_source(...) / d^2`.

Returns an area (m^2). For a point-source Sun this reduces to
`(2/3) (albedo/pi) R^2 [sin(a) + (pi - a) cos(a)]`.

# Arguments
- `target`, `sun`, `observer`: States (any units).
- `radius`: Sphere radius (any length unit).
- `albedo`: Lambertian (diffuse) reflectance, 0-1.
"""
function lambertian_extended_source(target, sun, observer, radius, albedo)
    r = to_unit(u"m", radius)
    rho = to_fraction(albedo)
    phase = solar_phase_angle(target, sun, observer)
    theta_s = solar_apparent_radius(target, sun)

    return (1 / 6) * (rho / pi) * r^2 * cos(theta_s) * (5 - cos(theta_s)) *
           (sin(phase) + (pi - phase) * cos(phase))
end
