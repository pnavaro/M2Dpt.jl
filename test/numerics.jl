# Numerics
const nx     = 94      :: Int64    # number of cells x
const ny     = nx      :: Int64    # number of cells y
const nt     = 54      :: Int64    # number time steps
const nout   = 100     :: Int64    # check residual each nout iteration
const noutp  = 2       :: Int64    # display graphics each nout time step
const niter  = 100000  :: Int64    # max nonlinear iterations
const epsi   = 1e-5    :: Float64  # non-linear tolerance
const tetp   = 0.5     :: Float64  # reduction of PT steps for pressure
const tetv   = 0.5     :: Float64  # reduction of PT steps for velocity
const tetT   = 0.5     :: Float64  # reduction of physical time step for temperature
const rel    = 1e-1    :: Float64  # relaxation factor for non-linear viscosity
const Vdamp  = 4.0     :: Float64  # velocity damping for momentum equations
const eta_b  = 1.0     :: Float64  # numerical compressibility
