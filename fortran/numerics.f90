module numerics 

    integer, parameter   :: nx     = 94      ! number of cells x
    integer, parameter   :: ny     = nx      ! number of cells y
    integer, parameter   :: nt     = 54      ! number time steps
    integer, parameter   :: nout   = 100     ! check residual each nout iteration
    integer, parameter   :: noutp  = 2       ! display graphics each nout time step
    integer, parameter   :: niter  = 1e5     ! max nonlinear iterations
    real(8), parameter   :: epsi   = 1e-5    ! non-linear tolerance
    real(8), parameter   :: tetp   = 0.5     ! reduction of PT steps for pressure
    real(8), parameter   :: tetv   = 0.5     ! reduction of PT steps for velocity
    real(8), parameter   :: tetT   = 0.5     ! reduction of physical time step for temperature
    real(8), parameter   :: rel    = 1e-1    ! relaxation factor for non-linear viscosity
    real(8), parameter   :: Vdamp  = 4.0     ! velocity damping for momentum equations
    real(8), parameter   :: eta_b  = 1.0     ! numerical compressibility

end module numerics
