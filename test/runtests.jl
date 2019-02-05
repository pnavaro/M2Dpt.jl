using Printf
using LinearAlgebra

include("../src/mesh.jl")
include("../src/fields.jl")
include("physics.jl")
include("numerics.jl")

# Pre-processing
const dampx = 1*(1-Vdamp/nx) # velocity damping for x-momentum equation
const dampy = 1*(1-Vdamp/ny) # velocity damping for y-momentum equation
const mpow  = -(1-1/n)/2     # exponent for strain rate dependent viscosity

include("vectorized_solve.jl")

mesh = Mesh( Lx, nx, Ly, ny)

fields = Fields( mesh, Vbc, r, Tamp )

@time vectorized_solve( mesh, fields )
