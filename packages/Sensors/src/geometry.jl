#=
Observation geometry built on SPICE. All inputs accept Unitful quantities;
all states are returned as `CartesianState`.

Requires the leapseconds, planetary ephemeris and PCK kernels to be loaded
(e.g. `Astrodynamics.furnsh_kernels()`).
=#

"""
    ephemeris_time(t)

Convert `t` to TDB seconds past J2000. Accepts a number (already ET seconds)
or a SPICE-parsable time string such as `"2026-10-10T23:00:00 UTC"`.
"""
ephemeris_time(t::Real) = Float64(t)
ephemeris_time(t::AbstractString) = str2et(t)

"""
    spice_state(target, et; frame="J2000", observer="EARTH", abcorr="NONE")

State of a SPICE body relative to `observer`, as a `CartesianState`.
"""
function spice_state(
    target::AbstractString,
    et;
    frame::AbstractString = "J2000",
    observer::AbstractString = "EARTH",
    abcorr::AbstractString = "NONE"
)
    x, _ = spkezr(target, ephemeris_time(et), frame, abcorr, observer)
    return CartesianState(x, u"km", u"km/s")
end

"""
    _body_fixed_to_inertial(r_bf, v_bf, et, frame)

Rotate an Earth-fixed (IAU_EARTH) position/velocity to `frame`, including the
transport term from Earth rotation.
"""
function _body_fixed_to_inertial(r_bf_km, v_bf_kms, et, frame)
    xform = sxform("IAU_EARTH", frame, ephemeris_time(et))
    x = xform * vcat(collect(r_bf_km), collect(v_bf_kms))
    return CartesianState(x, u"km", u"km/s")
end

"""
    GroundSite(; longitude, latitude, altitude)

Geodetic ground site on the Earth ellipsoid from the loaded PCK.

# Fields
- `longitude::Radians`: East longitude.
- `latitude::Radians`: Geodetic latitude.
- `altitude::Meters`: Height above the ellipsoid.
"""
struct GroundSite
    longitude::Radians
    latitude::Radians
    altitude::Meters
end

function GroundSite(; longitude, latitude, altitude)
    return GroundSite(
        to_unit(u"rad", longitude),
        to_unit(u"rad", latitude),
        to_unit(u"m", altitude)
    )
end

"""
    GroundSite(longitude, rho_cos_phi, rho_sin_phi)

Build a site from Minor Planet Center parallax constants (in Earth radii,
R = 6378.137 km), as listed at
https://minorplanetcenter.net/iau/lists/ObsCodes.html.
"""
function GroundSite(longitude::Quantity, rho_cos_phi::Real, rho_sin_phi::Real)
    lon = ustrip(u"rad", to_unit(u"rad", longitude))
    r_eq = 6378.137 * rho_cos_phi
    r_bf = [r_eq * cos(lon), r_eq * sin(lon), 6378.137 * rho_sin_phi]
    re, _, rp = bodvrd("EARTH", "RADII")
    lon_g, lat_g, alt_km = recgeo(r_bf, re, (re - rp) / re)
    return GroundSite(lon_g * u"rad", lat_g * u"rad", alt_km * u"km")
end

"""
    site_state(site, et; frame="J2000")

Inertial state of a ground site (Earth-centred).
"""
function site_state(site::GroundSite, et; frame::AbstractString = "J2000")
    re, _, rp = bodvrd("EARTH", "RADII")
    r_bf = georec(
        ustrip(u"rad", site.longitude),
        ustrip(u"rad", site.latitude),
        ustrip(u"km", site.altitude),
        re,
        (re - rp) / re
    )
    return _body_fixed_to_inertial(r_bf, zeros(3), et, frame)
end

"""
    site_zenith(site, et; frame="J2000")

Unit vector of the local (geodetic) vertical of `site` in `frame`.
"""
function site_zenith(site::GroundSite, et; frame::AbstractString = "J2000")
    lon = ustrip(u"rad", site.longitude)
    lat = ustrip(u"rad", site.latitude)
    up_bf = [cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)]
    rot = pxform("IAU_EARTH", frame, ephemeris_time(et))
    return SVector{3}(rot * up_bf)
end

"""
    geostationary_state(; longitude, et, frame="J2000")

State of an ideal geostationary object (fixed over `longitude` in the
IAU_EARTH frame, at radius `GEO_RADIUS`).
"""
function geostationary_state(; longitude, et, frame::AbstractString = "J2000")
    lon = ustrip(u"rad", to_unit(u"rad", longitude))
    a = ustrip(u"km", GEO_RADIUS)
    r_bf = [a * cos(lon), a * sin(lon), 0.0]
    return _body_fixed_to_inertial(r_bf, zeros(3), et, frame)
end

"""
    line_of_sight(target, observer)

Unit vector from `observer` to `target`.
"""
function line_of_sight(target, observer)
    rho = CartesianState(target).r - CartesianState(observer).r
    return SVector{3, Float64}(rho / norm(rho))
end

"""
    radec(target, observer)

Right ascension and declination of `target` seen from `observer`, in the
frame of the input states. Returns a tuple of `Radians`.
"""
function radec(target, observer)
    u = line_of_sight(target, observer)
    return (atan(u[2], u[1]) * u"rad", asin(clamp(u[3], -1.0, 1.0)) * u"rad")
end
