program m2dpt

    implicit none

    integer :: i
    integer :: j

    integer :: it
    integer :: iter
    real(8) :: time
    real(8) :: etac_phys
    real(8) :: err_fu
    real(8) :: err_fp
    real(8) :: err_ft
    real(8) :: dvxdx
    real(8) :: dvydy
    real(8) :: dtt
    real(8) :: dtauT
    real(8) :: dx 
    real(8), allocatable :: xc(:)
    real(8), allocatable :: xn(:)
    real(8) :: dy 
    real(8), allocatable :: yc(:) 
    real(8), allocatable :: yn(:) 

    real(8), allocatable :: Vx(:,:)
    real(8), allocatable :: Vy(:,:)
    real(8), allocatable :: T(:,:)
    real(8), allocatable :: P(:,:)
    real(8), allocatable :: etac(:,:)
    real(8), allocatable :: qx(:,:)
    real(8), allocatable :: qy(:,:)
    real(8), allocatable :: dVxdtauVx(:,:)
    real(8), allocatable :: dVydtauVy(:,:)
    real(8), allocatable :: dVxdtauVx0(:,:)
    real(8), allocatable :: dVydtauVy0(:,:)

    integer, parameter   :: n      = 3          ! stress exponent for power law rheology
    real(8), parameter   :: Vbc    = 66.4437    ! boundary velocity difference
    real(8), parameter   :: Lx     = 0.86038    ! domain size x
    real(8), parameter   :: Ly     = Lx         ! domain size y
    real(8), parameter   :: T0     = 49.3269/n  ! initial temperature
    real(8), parameter   :: r      = 0.0737     ! radius of initial perturbation
    real(8), parameter   :: Tamp   = 0.1*T0     ! amplitude of initial perturbation

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

    real(8), parameter   :: dampx = 1*(1-Vdamp/nx)   ! velocity damping for x-momentum equation
    real(8), parameter   :: dampy = 1*(1-Vdamp/ny)   ! velocity damping for y-momentum equation
    real(8), parameter   :: mpow  = -(1d0-1d0/n)/2d0 ! exponent for strain rate dependent viscosity

    real(8), allocatable :: Vx_exp  (:,:)
    real(8), allocatable :: Vy_exp  (:,:)
    real(8), allocatable :: etav    (:,:) 
    real(8), allocatable :: Exyv    (:,:) 
    real(8), allocatable :: Exyc    (:,:)
    real(8), allocatable :: divV    (:,:)
    real(8), allocatable :: Exxc    (:,:)
    real(8), allocatable :: Eyyc    (:,:)
    real(8), allocatable :: Eii2    (:,:)
    real(8), allocatable :: Sxx     (:,:)
    real(8), allocatable :: Syy     (:,:)
    real(8), allocatable :: Txy     (:,:)
    real(8), allocatable :: Hs      (:,:)
    real(8), allocatable :: To      (:,:)
    real(8), allocatable :: dPdtauP (:,:)
    real(8), allocatable :: dTdtauT (:,:)
    real(8), allocatable :: dtauP   (:,:)
    real(8), allocatable :: dtauVx  (:,:)
    real(8), allocatable :: dtauVy  (:,:)


    dx = Lx/nx
    dy = Ly/ny

    allocate(xn(nx+1))
    allocate(xc(nx))
    allocate(yn(ny+1))
    allocate(yc(ny))

    do i = 1, nx+1
       xn(i) = (i-1)*dx
    end do
    do i = 1, nx
       xc(i) = (i-0.5)*dx
    end do
    do j = 1, ny+1
       yn(j) = (j-1)*dy
    end do
    do j = 1, ny
       yc(j) = (j-0.5)*dy
    end do

    allocate(Vx(nx+1,ny))
    allocate(Vy(nx,ny+1))

    do j = 1, ny+1
        do i = 1, nx+1
            if (j <= ny) Vx(i,j) =  Vbc * xn(i) /Lx
            if (i <= nx) Vy(i,j) = -Vbc * yn(j) /Ly
        end do
    end do

    allocate(T         (nx  ,ny  )); T    = 0d0
    allocate(P         (nx  ,ny  )); P    = 0d0
    allocate(etac      (nx  ,ny  )); etac = 1d0
    allocate(qx        (nx+1,ny  )); qx   = 0d0
    allocate(qy        (nx  ,ny+1)); qy   = 0d0
    allocate(dVxdtauVx (nx-1,ny  )); dVxdtauVx  = 0d0
    allocate(dVydtauVy (nx  ,ny-1)); dVydtauVy  = 0d0
    allocate(dVxdtauVx0(nx-1,ny  )); dVxdtauVx0 = 0d0
    allocate(dVydtauVy0(nx  ,ny-1)); dVydtauVy0 = 0d0

    do j = 1, ny
        do i = 1, nx
            if (xc(i)**2 + yc(j)**2  < r**2) then
                T(i,j) = Tamp        ! initial temperature pertubation
            else
                T(i,j) = 0d0
            end if
        end do
    end do

    time   = 0.0

    dtT = tetT*1/4.1*min(dx,dy)**2 ! explicit timestep for 2D diffusion
    
    allocate(Vx_exp  (nx+1,ny+2)); Vx_exp  = 0d0
    allocate(Vy_exp  (nx+2,ny+1)); Vy_exp  = 0d0
    allocate(etav    (nx+1,ny+1)); etav    = 0d0 
    allocate(Exyv    (nx+1,ny+1)); Exyv    = 0d0 
    allocate(Exyc    (nx,ny))    ; Exyc    = 0d0 
    allocate(divV    (nx,ny))    ; divV    = 0d0
    allocate(Exxc    (nx,ny))    ; Exxc    = 0d0
    allocate(Eyyc    (nx,ny))    ; Eyyc    = 0d0
    allocate(Eii2    (nx,ny))    ; Eii2    = 0d0
    allocate(Sxx     (nx,ny))    ; Sxx     = 0d0
    allocate(Syy     (nx,ny))    ; Syy     = 0d0
    allocate(Txy     (nx+1,ny+1)); Txy     = 0d0
    allocate(Hs      (nx,ny))    ; Hs      = 0d0
    allocate(To      (nx,ny))    ; To      = 0d0
    allocate(dPdtauP (nx,ny))    ; dPdtauP = 0d0
    allocate(dTdtauT (nx,ny))    ; dTdtauT = 0d0
    allocate(dtauP   (nx,ny))    ; dtauP   = 0d0
    allocate(dtauVx  (nx-1,ny))  ; dtauVx  = 0d0
    allocate(dtauVy  (nx,ny-1))  ; dtauVy  = 0d0

    do it = 1,nt ! Physical timesteps

        To = T ! temperature from previous step (for backward-Euler integration)

        time = time + dtT ! update physical time

        do iter = 1,niter ! Pseudo-Transient cycles

            ! used for damping x momentum residuals
            dVxdtauVx0 = dVxdtauVx + dampx * dVxdtauVx0  

            ! used for damping y momentum residuals
            dVydtauVy0 = dVydtauVy + dampy * dVydtauVy0 

            !  Kinematics

            do i = 1,nx+1
                Vx_exp(i,1)  = Vx(i,1)
                do j = 1,ny
                    Vx_exp(i,j+1) = Vx(i,j)
                end do
                Vx_exp(i,ny+2) = Vx(i,ny)
            end do
                
            do j = 1,ny+1
                Vy_exp(1,j) = Vy(1,j)
                do i = 1,nx
                    Vy_exp(i+1,j) = Vy(i,j)
                end do
                Vy_exp(nx+2,j) = Vy(nx,j)
            end do

            do j=1,ny
            do i=1,nx

                dVxdx = (Vx(i+1,j)-Vx(i,j))/dx
                dVydy = (Vy(i,j+1)-Vy(i,j))/dy

                divV(i,j) = dVxdx + dVydy

                Exxc(i,j) = dVxdx - 0.5d0 * (dVxdx + dVydy)

                Eyyc(i,j) = dVydy - 0.5d0 * (dVxdx + dVydy)

            end do
            end do



            do j=1,ny+1
            do i=1,nx+1
                Exyv(i,j) = 0.5d0*(  (Vx_exp(i,j+1)-Vx_exp(i,j))/dy  &
                                   + (Vy_exp(i+1,j)-Vy_exp(i,j))/dx )
            end do
            end do

            do j=1,ny
            do i=1,nx
                Exyc(i,j) = 0.25d0*(Exyv(i,j  )+Exyv(i+1,j  ) &
                                   +Exyv(i,j+1)+Exyv(i+1,j+1))
            end do
            end do

            do j=1,ny
            do i=1,nx
                Eii2(i,j) = 0.5d0*(Exxc(i,j)**2 + Eyyc(i,j)**2) + Exyc(i,j)**2 ! strain rate invariant
            end do 
            end do 


            ! ------ Rheology

            do j=1,ny
            do i=1,nx

                 ! physical viscosity
                 etac_phys = Eii2(i,j)**mpow*exp( -T(i,j)*(1 / (1 + T(i,j)/T0)) ) 

                 ! numerical shear viscosity
                 etac(i,j) = exp(rel*log(etac_phys) + (1-rel)*log(etac(i,j)))

            end do
            end do

            ! expand viscosity fom cell centroids to vertices
            etav = 0.0

            do j = 2,ny
            do i = 2,nx
                etav(i,j) = 0.25*(etac(i-1,j-1)+etac(i,j-1) &
                                 +etac(i-1,j  )+etac(i,j  ))
            end do
            end do

            do i = 1,nx+1
                etav(i,1   ) = etav(i, 2)
                etav(i,ny+1) = etav(i,ny)
            end do
            do j = 1,ny+1
                etav(1,   j) = etav(2, j)
                etav(nx+1,j) = etav(nx,j)
            end do

            ! ------ Pseudo-Time steps ------

            dtauP = tetp *  4.1 / min(nx,ny)*etac*(1.0+eta_b)

            do j=1,ny
            do i=1,nx-1
                dtauVx(i,j) = tetv * 1/4.1 * (min(dx,dy)**2 / ( & 
                          0.5*(etac(i+1,j) + etac(i,j))))/(1+eta_b)
            end do
            end do

            do j=1,ny-1
            do i=1,nx
                dtauVy(i,j) = tetv * 1/4.1 * (min(dx,dy)**2 / ( &
                           0.5*(etac(i,j+1) + etac(i,j))))/(1+eta_b)
            end do
            end do

            dtauT = tetT * 1/4.1 * min(dx,dy)**2

            ! ------ Fluxes

            do j=1,ny
            do i=2,nx
                qx(i,j) = -(T(i,j)-T(i-1,j))/dx
            end do
            end do

            do j=2,ny
            do i=1,nx
                qy(i,j) = -(T(i,j)-T(i,j-1))/dy
            end do
            end do

            do j=1,ny
            do i=1,nx
                Sxx(i,j) = -P(i,j) + 2 * etac(i,j) * (Exxc(i,j) + eta_b*divV(i,j))
                Syy(i,j) = -P(i,j) + 2 * etac(i,j) * (Eyyc(i,j) + eta_b*divV(i,j))
            end do
            end do

            do j=1,ny+1
            do i=1,nx+1
                Txy(i,j) = 2 * etav(i,j)  * Exyv(i,j)
            end do
            end do

            do j=1,ny
            do i=1,nx
                Hs(i,j) = 4 * etac(i,j) * Eii2(i,j)
            end do
            end do

            ! ------ Residuals

            do j=1,ny
            do i=1,nx-1
                dVxdtauVx(i,j) = (Txy(i+1,j+1)-Txy(i+1,j))/dy + (Sxx(i+1,j)-Sxx(i,j))/dx
            end do
            end do

            do j=1,ny-1
            do i=1,nx
                dVydtauVy(i,j) = (Txy(i+1,j+1)-Txy(i,j+1))/dx + (Syy(i,j+1)-Syy(i,j))/dy
            end do
            end do

             
            do j=1,ny
            do i=1,nx
                dPdtauP(i,j) = - divV(i,j)
                dTdtauT(i,j) = ((To(i,j)-T(i,j))/dtT  &
                               - ((qx(i+1,j)-qx(i,j))/dx + (qy(i,j+1)-qy(i,j))/dy)  &
                               + Hs(i,j))
            end do
            end do
            ! ------ Updates

            ! update with damping
            Vx(2:nx,:) = Vx(2:nx,:) + dtauVx * (dVxdtauVx + dampx*dVxdtauVx0) 
            Vy(:,2:ny) = Vy(:,2:ny) + dtauVy * (dVydtauVy + dampy*dVydtauVy0)
            P          = P + dtauP  * dPdtauP
            T          = T + dtauT  * dTdtauT

            if (mod(iter,nout) ==0) then! Check

                err_fu = sqrt(sum(dVxdtauVx**2)+sum(dVydtauVy**2))/(nx*(ny-1)+ny*(ny-1))
                err_fp = sqrt(sum(dPdtauP**2))/(nx*ny)
                err_fT = sqrt(sum(dTdtauT**2))/(nx*ny)

                if (max(err_fu, err_fp, err_fT) < epsi) then
                    exit
                end if

                print*,"-------------------------"
                print" ('iter  = ',i5   )", iter  
                print" ('f_{u} = ',e10.3 )", err_fu
                print" ('f_{p} = ',e10.3 )", err_fp
                print" ('f_{T} = ',e10.3 )", err_fT

            end if

        end do

    end do

end program m2dpt
