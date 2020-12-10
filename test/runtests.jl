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

include("simd_test.jl")
include("vectorized_solve.jl")
include("vectorized_solve_Diff.jl")
include("solve_with_loops.jl")
include("solve_with_loops_LoopVec.jl")

smalln = [4, 40][1] # 36 also works but slowly, numbers 4<smalln<36 produce NaNs(?!)
println("vectorized_solve")
mesh = Mesh( Lx, nx, Ly, ny)
fields = Fields( mesh, Vbc, r, Tamp )
mesh_small = Mesh( Lx, smalln, Ly, smalln)
fields_small = Fields( mesh_small, Vbc, r, Tamp )
@time vectorized_solve( mesh_small, fields_small)
f1s = deepcopy(fields_small)
@time vectorized_solve( mesh, fields )
f1 = deepcopy(fields)

println("vectorized_solve with Diff-trick")
mesh = Mesh( Lx, nx, Ly, ny)
fields = Fields( mesh, Vbc, r, Tamp )
mesh_small = Mesh( Lx, smalln, Ly, smalln)
fields_small = Fields( mesh_small, Vbc, r, Tamp )
@time vectorized_solve_v2( mesh_small, fields_small)
f2s = deepcopy(fields_small)
@time vectorized_solve_v2( mesh, fields )
f2 = deepcopy(fields)

println("loop")
mesh = Mesh( Lx, nx, Ly, ny)
fields = Fields( mesh, Vbc, r, Tamp )
mesh_small = Mesh( Lx, smalln, Ly, smalln)
fields_small = Fields( mesh_small, Vbc, r, Tamp )
@time solve_with_loops( mesh_small, fields_small)
f3s = deepcopy(fields_small)
@time solve_with_loops( mesh, fields )
f3 = deepcopy(fields)

println("loop with LoopVectorizations.jl")
mesh = Mesh( Lx, nx, Ly, ny)
fields = Fields( mesh, Vbc, r, Tamp )
mesh_small = Mesh( Lx, smalln, Ly, smalln)
fields_small = Fields( mesh_small, Vbc, r, Tamp )
@time solve_with_loops_v2( mesh_small, fields_small)
f4s = deepcopy(fields_small)
@time solve_with_loops_v2( mesh, fields )
f4 = deepcopy(fields)


# check we get the same
@assert(all(f1s.T.==f2s.T))
@assert(all(f1.T.==f2.T))
@assert(all(f1.T.≈f3.T))
@assert(all(f1.T.≈f4.T))
;
