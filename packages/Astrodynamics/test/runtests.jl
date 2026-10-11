using Astrodynamics
using Test
using StaticArrays
using LinearAlgebra

@testset "Coordinate Conversion Tests" begin
    # Define a sample observer state [x, y, z, vx, vy, vz]
    observer = SVector(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
    
    # Define a sample spacecraft state in Cartesian coordinates
    state_cart = SVector(7100.0, 100.0, 50.0, -1.0, 7.4, 0.1)

    @testset "cart2sph output format" begin
        state_sph = cart2sph(state_cart, observer)
        
        # Ensure output is an SVector of length 6
        @test state_sph isa SVector{6, Float64}
        @test length(state_sph) == 6
        
        # Range must be positive
        @test state_sph[1] > 0.0
    end

    @testset "Round-trip Conversion (Cartesian -> Spherical -> Cartesian)" begin
        # Convert Cartesian to Spherical
        state_sph = cart2sph(state_cart, observer)
        
        # Convert back to Cartesian
        state_cart_recovered = sph2cart(state_sph, observer)
        
        # Check that we recover the original state within floating-point tolerance
        @test state_cart_recovered ≈ state_cart atol=1e-10
    end

    @testset "Round-trip Conversion (Spherical -> Cartesian -> Spherical)" begin
        # Define a sample state in spherical coordinates: [r, ra, dec, r_dot, ra_dot, dec_dot]
        state_sph_original = SVector(1000.0, 1.2, 0.3, 0.5, 0.01, -0.02)
        
        # Convert Spherical to Cartesian
        state_cart_interp = sph2cart(state_sph_original, observer)
        
        # Convert back to Spherical
        state_sph_recovered = cart2sph(state_cart_interp, observer)
        
        # Check that we recover the original spherical state
        @test state_sph_recovered ≈ state_sph_original atol=1e-10
    end
end