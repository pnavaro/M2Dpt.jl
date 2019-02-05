module physics

    integer, parameter   :: n      = 3          ! stress exponent for power law rheology
    real(8), parameter   :: Vbc    = 66.4437    ! boundary velocity difference
    real(8), parameter   :: Lx     = 0.86038    ! domain size x
    real(8), parameter   :: Ly     = Lx         ! domain size y
    real(8), parameter   :: T0     = 49.3269/n  ! initial temperature
    real(8), parameter   :: r      = 0.0737     ! radius of initial perturbation
    real(8), parameter   :: Tamp   = 0.1*T0     ! amplitude of initial perturbation

end module physics
