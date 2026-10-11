"""
    propagate(dynamics, initial_state, tspan; params=nothing, abstol=1e-9, reltol=1e-9, solver=Vern9(), callback=nothing, kwargs...)

Propagation wrapper for any dynamic function using OrdinaryDiffEq.jl.

# Arguments
- `dynamics::Function`: The dynamics function to propagate.
- `initial_state::SVector`: The initial state vector.
- `tspan::Tuple`: The time span for the propagation.
- `params`: Optional parameters to pass to the dynamics function.
- `abstol`: Absolute tolerance for the ODE solver (default is 1e-9).
- `reltol`: Relative tolerance for the ODE solver (default is 1e-9).
- `solver`: The ODE solver to use (default is Vern9()).
- `callback`: Optional callback function for the ODE solver.
- `kwargs...`: Additional keyword arguments to pass to the ODE solver.

# Returns
- `sol`: The solution object returned by the ODE solver.
"""
function propagate(
    dynamics::Function,
    initial_state::SVector,
    tspan::Tuple;
    params = nothing,
    abstol = 1e-9,
    reltol = 1e-9,
    solver = Vern9(),
    callback = nothing,
    kwargs...
)
    prob = ODEProblem(dynamics, initial_state, tspan, params)
    sol = solve(
        prob, solver,
        abstol = abstol, reltol = reltol,
        callback = callback, kwargs...
    )
    return sol
end

"""
    CelestialBody{F}
A struct representing a celestial body with its gravitational parameter and a spline function for its state.

# Fields
- `gm::Float64`: The gravitational parameter of the celestial body.
- `spline::F`: A function that returns a 6-element position/velocity vector at a given epoch `t`.
"""
struct CelestialBody{F}
    gm::Float64
    spline::F  # Function returning 6-element position/velocity vector at epoch t
end

"""
    NBodyParams{T}
A struct for holding simulation parameters for an N-body orbit mechanics problem.

# Fields
- `initial_et::Float64`: The initial ephemeris time.
- `central_gm::Float64`: The gravitational parameter of the central body.
"""
struct NBodyParams{T}
    initial_et::Float64
    central_gm::Float64
    bodies::T  # Tuple or Vector of CelestialBody structs
end

"""
    n_body_spice(state, params::NBodyParams, t)
Compute the time derivative of the state vector for an N-body problem using SPICE data.

# Arguments
- `state::SVector`: The current state vector (position and velocity).
- `params::NBodyParams`: The parameters for the N-body simulation.
- `t::Float64`: The current time in seconds since the initial epoch.
"""
function n_body_spice(state, params::NBodyParams, t)
    # Extract 3D position and velocity directly into StaticArrays
    r_sc = SVector(state[1], state[2], state[3])
    v_sc = SVector(state[4], state[5], state[6])
    
    t_et = t + params.initial_et

    # 1. Primary body central gravity acceleration
    r_sc_mag = norm(r_sc)
    a_sc = -params.central_gm * r_sc / r_sc_mag^3

    # 2. Add third-body perturbations from N bodies
    for body in params.bodies
        body_state = body.spline(t_et)
        r_body = SVector(body_state[1], body_state[2], body_state[3])
        
        r_body_sc = r_body - r_sc
        
        r_body_sc_mag = norm(r_body_sc)
        r_body_mag = norm(r_body)

        # Indirect + Direct third-body gravity formulation
        a_sc += body.gm * (r_body_sc / r_body_sc_mag^3 - r_body / r_body_mag^3)
    end

    # Explicit 6-element SVector construction
    return SVector(v_sc[1], v_sc[2], v_sc[3], a_sc[1], a_sc[2], a_sc[3])
end

