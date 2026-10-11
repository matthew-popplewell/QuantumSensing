"""
    cart2sph(state::SVector{6, Float64}, observer::SVector{6, Float64})
Convert a state vector from Cartesian coordinates to spherical coordinates relative to an observer.

# Arguments
- `state::SVector{6, Float64}`: The state vector in Cartesian coordinates.
- `observer::SVector{6, Float64}`: The observer's state vector in Cartesian coordinates.
"""
function cart2sph(state::SVector{6, Float64}, observer::SVector{6, Float64})
    # Extract position and velocity vectors
    r_sc = SVector(state[1], state[2], state[3])
    v_sc = SVector(state[4], state[5], state[6])
    r_obs = SVector(observer[1], observer[2], observer[3])
    v_obs = SVector(observer[4], observer[5], observer[6])

    # Compute relative position and velocity
    r_rel = r_sc - r_obs
    v_rel = v_sc - v_obs

    # Compute spherical coordinates
    r_mag = norm(r_rel)
    los = r_rel / r_mag  # Line-of-sight unit vector
    ra = atan(los[2], los[1])
    dec = asin(los[3])
    
    range_rate = dot(r_rel, v_rel) / r_mag
    los_rate = (v_rel - range_rate * los) / r_mag
    
    sin_ra, cos_ra = sin(ra), cos(ra)
    sin_dec, cos_dec = sin(dec), cos(dec)
    dec_rate = los_rate[3] / cos_dec
    ra_rate = (los_rate[2] + dec_rate * sin_ra * sin_dec) / (cos_ra * cos_dec)

    return SVector(r_mag, ra, dec, range_rate, ra_rate, dec_rate)
end

"""
    sph2cart(state::SVector{6, Float64}, observer::SVector{6, Float64})
Convert a state vector from spherical coordinates to Cartesian coordinates relative to an observer.

# Arguments
- `state::SVector{6, Float64}`: The state vector in spherical coordinates.
- `observer::SVector{6, Float64}`: The observer's state vector in Cartesian coordinates.
"""
function sph2cart(state::SVector{6, Float64}, observer::SVector{6, Float64})
    # Extract spherical coordinates and rates
    r_mag = state[1]
    ra = state[2]
    dec = state[3]
    range_rate = state[4]
    ra_rate = state[5]
    dec_rate = state[6]

    sin_ra, cos_ra = sin(ra), cos(ra)
    sin_dec, cos_dec = sin(dec), cos(dec)

    # Compute line-of-sight unit vector
    los = SVector(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    los_rate = SVector(
        -ra_rate * sin_ra * cos_dec - dec_rate * cos_ra * sin_dec,
        ra_rate * cos_ra * cos_dec - dec_rate * sin_ra * sin_dec,
        dec_rate * cos_dec
    )

    # Compute relative position and velocity
    r_rel = r_mag * los
    v_rel = range_rate * los + r_mag * los_rate

    # Compute absolute position and velocity
    r_obs = SVector(observer[1], observer[2], observer[3])
    v_obs = SVector(observer[4], observer[5], observer[6])
    
    r_sc = r_obs + r_rel
    v_sc = v_obs + v_rel

    return SVector(r_sc[1], r_sc[2], r_sc[3], v_sc[1], v_sc[2], v_sc[3])
end