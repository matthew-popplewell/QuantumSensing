#=
Unit handling.

Every public constructor and function accepts any Unitful quantity with the
right dimension (km, m, AU, μm, arcsec, deg, ...) and converts it once to the
canonical SI unit stored below. Internally everything is concretely typed, so
Unitful adds no runtime cost and still checks the dimensional consistency of
the radiometry.
=#

# Canonical, concrete storage types.
const Meters = typeof(1.0u"m")
const MetersPerSecond = typeof(1.0u"m/s")
const SquareMeters = typeof(1.0u"m^2")
const Seconds = typeof(1.0u"s")
const PerSecond = typeof(1.0u"s^-1")
const Radians = typeof(1.0u"rad")
const SpectralFlux = typeof(1.0u"W/m^3")  # F_lambda, W m^-2 m^-1

"""
    to_unit(unit, x)

Convert `x` to `unit` as a Float64 quantity.

Throws a `DimensionError` if `x` has incompatible units, and an
`ArgumentError` if `x` is a plain number but `unit` is not dimensionless.
Plain numbers are rejected for angles as well, because "0.5" could mean
radians, degrees or arcseconds.
"""
to_unit(unit::Unitful.Units, x::Quantity) = float(uconvert(unit, x))

function to_unit(unit::Unitful.Units, x::Real)
    throw(ArgumentError(
        "Got the plain number $x where a quantity in units compatible " *
        "with $unit was expected. Attach units, e.g. `$(x)u\"$unit\"`."
    ))
end

"""
    to_fraction(x)

Convert a dimensionless fraction to Float64. Accepts `0.6`, `60u"percent"`,
etc.
"""
to_fraction(x::Real) = Float64(x)
to_fraction(x::Quantity) = Float64(uconvert(NoUnits, x))

"""
    CartesianState(r, v)
    CartesianState(x)
    CartesianState(x, length_unit, velocity_unit)

Cartesian position and velocity, stored in m and m/s.

# Constructors
- `CartesianState(r, v)`: `r` and `v` are 3-vectors of any length / velocity
  units (they may be mixed, e.g. `[1u"AU", 0u"km", 0u"m"]`).
- `CartesianState(x)`: `x` is a 6-vector of quantities `[r; v]`.
- `CartesianState(x, u"km", u"km/s")`: `x` is a plain 6-vector (e.g. SPICE
  output) whose units you state explicitly.

# Fields
- `r::SVector{3, Meters}`: Position.
- `v::SVector{3, MetersPerSecond}`: Velocity.
"""
struct CartesianState
    r::SVector{3, Meters}
    v::SVector{3, MetersPerSecond}
end

function CartesianState(r::AbstractVector, v::AbstractVector)
    length(r) == 3 && length(v) == 3 ||
        throw(ArgumentError("Position and velocity must have 3 components."))
    return CartesianState(
        SVector{3, Meters}(to_unit.(u"m", r)),
        SVector{3, MetersPerSecond}(to_unit.(u"m/s", v))
    )
end

function CartesianState(x::AbstractVector)
    length(x) == 6 ||
        throw(ArgumentError("A state vector must have 6 components."))
    return CartesianState(x[1:3], x[4:6])
end

function CartesianState(
    x::AbstractVector{<:Real},
    length_unit::Unitful.Units,
    velocity_unit::Unitful.Units
)
    length(x) == 6 ||
        throw(ArgumentError("A state vector must have 6 components."))
    return CartesianState(x[1:3] .* length_unit, x[4:6] .* velocity_unit)
end

CartesianState(x::CartesianState) = x

function CartesianState(x::AbstractVector{<:Real})
    throw(ArgumentError(
        "Got a plain-number state vector. State its units explicitly, e.g. " *
        "`CartesianState(x, u\"km\", u\"km/s\")` for SPICE output."
    ))
end

"""
    state_vector(s; length_unit=u"km", velocity_unit=u"km/s")

Strip units from a `CartesianState` into a plain `SVector{6, Float64}`, e.g.
for use with `Astrodynamics.propagate` (km, km/s convention).
"""
function state_vector(
    s::CartesianState;
    length_unit::Unitful.Units = u"km",
    velocity_unit::Unitful.Units = u"km/s"
)
    return SVector{6, Float64}(
        ustrip.(length_unit, s.r)..., ustrip.(velocity_unit, s.v)...
    )
end

"""
    drift(s, dt)

Straight-line propagation `r + v dt`. Adequate over a single exposure: for a
GEO target the neglected gravitational term is ~11 m over 10 s (~0.06 arcsec
at GEO range).
"""
function drift(s::CartesianState, dt)
    dt_s = to_unit(u"s", dt)
    return CartesianState(s.r + s.v * dt_s, s.v)
end

function Base.show(io::IO, s::CartesianState)
    r = round.(ustrip.(u"km", s.r), digits = 3)
    v = round.(ustrip.(u"km/s", s.v), digits = 6)
    print(io, "CartesianState(r = $r km, v = $v km/s)")
end
