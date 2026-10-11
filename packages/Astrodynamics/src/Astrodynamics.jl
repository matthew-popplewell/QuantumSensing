module Astrodynamics

using SPICE
using StaticArrays
using OrdinaryDiffEq
using LinearAlgebra

include("dynamics.jl")
include("spice.jl")
include("transformations.jl")

# dynamics.jl
export propagate
export CelestialBody
export NBodyParams
export n_body_spice

# spice.jl
export furnsh_kernels
export build_spline

# transformations.jl
export cart2sph
export sph2cart

end # module Astrodynamics
