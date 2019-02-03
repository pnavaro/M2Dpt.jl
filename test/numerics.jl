# Numerics
const nx     = 94                 # number of cells x
const ny     = nx                 # number of cells y
const nt     = 54                 # number time steps
const nout   = 100                # check residual each nout iteration
const noutp  = 2                  # display graphics each nout time step
const niter  = 1e5                # max nonlinear iterations
const epsi   = 1e-5               # non-linear tolerance
const tetp   = 0.5                # reduction of PT steps for pressure
const tetv   = 0.5                # reduction of PT steps for velocity
const tetT   = 0.5                # reduction of physical time step for temperature
const rel    = 1e-1               # relaxation factor for non-linear viscosity
const Vdamp  = 4.0                # velocity damping for momentum equations
const eta_b  = 1.0                # numerical compressibility
