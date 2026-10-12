#=
Visibility constraints for ground-based observations.
=#

"""
    elevation(site, et, target; frame="J2000")

Elevation of `target` above the local horizon of `site`.
"""
function elevation(site::GroundSite, et, target; frame::AbstractString = "J2000")
    u = line_of_sight(target, site_state(site, et; frame = frame))
    return asin(clamp(dot(u, site_zenith(site, et; frame = frame)), -1.0, 1.0)) *
           u"rad"
end

"""
    in_earth_shadow(target, sun)

Cylindrical Earth-shadow test. `target` and `sun` are Earth-centred states.
"""
function in_earth_shadow(target, sun)
    r = CartesianState(target).r
    s_hat = CartesianState(sun).r / norm(CartesianState(sun).r)
    along = dot(r, s_hat)
    along >= zero(along) && return false
    re = bodvrd("EARTH", "RADII")[1] * u"km"
    return norm(r - along * s_hat) < re
end

"""
    illumination_fraction(target, sun)

Fraction of the solar disk visible from `target` past the Earth (1 = full
Sun, 0 = umbra), from the overlapping-disk geometry in Robinson (2023)
Eq. 2.95-2.97 (after Krag 2003). `target` and `sun` are Earth-centred
states. Uses the Earth's equatorial radius.
"""
function illumination_fraction(target, sun)
    r_obj = CartesianState(target).r
    to_sun = CartesianState(sun).r - r_obj
    r_sun = bodvrd("SUN", "RADII")[1] * u"km"
    r_earth = bodvrd("EARTH", "RADII")[1] * u"km"

    eps_s = asin(uconvert(NoUnits, r_sun / norm(to_sun)))
    eps_e = asin(uconvert(NoUnits, r_earth / norm(r_obj)))
    eps = acos(clamp(dot(-r_obj / norm(r_obj), to_sun / norm(to_sun)), -1.0, 1.0))

    eps >= eps_s + eps_e && return 1.0
    eps <= eps_e - eps_s && return 0.0

    eps_1 = (eps^2 - eps_e^2 + eps_s^2) / (2 * eps)
    eps_2 = (eps^2 + eps_e^2 - eps_s^2) / (2 * eps)
    a_1 = pi * eps_s^2 - eps_s^2 * acos(clamp(eps_1 / eps_s, -1.0, 1.0)) +
          eps_1 * sqrt(max(eps_s^2 - eps_1^2, 0.0))
    a_2 = eps_e^2 * acos(clamp(eps_2 / eps_e, -1.0, 1.0)) -
          eps_2 * sqrt(max(eps_e^2 - eps_2^2, 0.0))
    return clamp((a_1 - a_2) / (pi * eps_s^2), 0.0, 1.0)
end
