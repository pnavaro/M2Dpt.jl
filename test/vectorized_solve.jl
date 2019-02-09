function vectorized_solve( m :: Mesh, f :: Fields )

    time   = 0.0
    dx, dy = mesh.dx, mesh.dy

    dtT = tetT*1/4.1*min(dx,dy)^2 # explicit timestep for 2D diffusion
    
    errs  = []

    etav    = zeros(Float64,(nx+1,ny+1)) 
    Exyv    = zeros(Float64,(nx+1,ny+1)) 
    Exyc    = zeros(Float64,(nx,ny)) 
    divV    = zeros(Float64,(nx,ny)) 
    Exxc    = zeros(Float64,(nx,ny)) 
    Eyyc    = zeros(Float64,(nx,ny)) 
    Eii2    = zeros(Float64,(nx,ny)) 
    Sxx     = zeros(Float64,(nx,ny)) 
    Syy     = zeros(Float64,(nx,ny)) 
    Txy     = similar(Exyv)
    Hs      = zeros(Float64,(nx,ny)) 
    To      = similar(f.T)
    dPdtauP = similar(divV)
    dTdtauT = similar(Hs)
    dtauP   = similar(f.etac)
    dtauVx  = zeros(Float64,(nx-1,ny))
    dtauVy  = zeros(Float64,(nx,ny-1))

    for it = 1:nt # Physical timesteps

        # temperature from previous step (for backward-Euler integration)
        To .= f.T 

        @show time += dtT # update physical time

        for iter = 1:niter # Pseudo-Transient cycles

            # used for damping x momentum residuals
            f.dVxdtauVx0 .= f.dVxdtauVx .+ dampx .* f.dVxdtauVx0  

            # used for damping y momentum residuals
            f.dVydtauVy0 .= f.dVydtauVy .+ dampy .* f.dVydtauVy0 

            #  Kinematics
            f.Vx[  :,  1] .= f.Vx[    :,    2]
            f.Vx[  :,end] .= f.Vx[    :,end-1]
            f.Vy[  1,  :] .= f.Vy[    2,    :]
            f.Vy[end,  :] .= f.Vy[end-1,    :]

            divV .= (diff(f.Vx[:,2:end-1],dims=1) ./ dx 
                  .+ diff(f.Vy[2:end-1,:],dims=2) ./ dy)

            Exxc .= diff(f.Vx[:,2:end-1],dims=1) ./ dx .- 1/2 .* divV
            Eyyc .= diff(f.Vy[2:end-1,:],dims=2) ./ dy .- 1/2 .* divV

            Exyv .= 0.5 .* (diff(f.Vx,dims=2)./dy .+ diff(f.Vy,dims=1)./dx)

  @views    Exyc .= 0.25 .* ( Exyv[1:end-1,1:end-1] .+ 
                              Exyv[2:end  ,1:end-1] .+ 
                              Exyv[1:end-1,2:end  ] .+ 
                              Exyv[2:end  ,2:end  ])

            # strain rate invariant
            Eii2 .= 0.5 .* (Exxc.^2 .+ Eyyc.^2) .+ Exyc.^2 

            # ------ Rheology
            # physical viscosity
            etac_phys = Eii2.^mpow .* exp.( -f.T.*(1 ./ (1 .+ f.T./T0)) ) 

            # numerical shear viscosity
            f.etac .= exp.(rel .* log.(etac_phys) .+ (1-rel) .* log.(f.etac)) 

            # expand viscosity fom cell centroids to vertices
            fill!(etav,0.0)

  @views    etav[2:end-1,2:end-1] .= 0.25*(f.etac[1:end-1,1:end-1] 
                                         + f.etac[2:end,2:end] 
                                         + f.etac[1:end-1,2:end] 
                                         + f.etac[2:end,1:end-1])

  @views    etav[:      ,[1 end]] .= etav[:        ,[2 end-1]]
  @views    etav[[1 end],:      ] .= etav[[2 end-1],:        ]

            # ------ Pseudo-Time steps ------

            dtauP   .= tetp *  4.1 / min(nx,ny)*f.etac*(1.0+eta_b)

  @views    dtauVx  .= tetv * 1/4.1 * (min(dx,dy)^2 ./ ( 
                      0.5.*(   f.etac[2:end,:] 
                           .+  f.etac[1:end-1,:]) ))./(1+eta_b)

  @views    dtauVy  .= tetv * 1/4.1 * (min(dx,dy)^2 ./ (
                       0.5.*(  f.etac[:,2:end] 
                           .+  f.etac[:,1:end-1]) ))./(1+eta_b)

            dtauT    = tetT * 1/4.1 * min(dx,dy)^2

            # ------ Fluxes

            f.qx[2:end-1,:] .= .- diff(f.T,dims=1) ./ dx
            f.qy[:,2:end-1] .= .- diff(f.T,dims=2) ./ dy

            Sxx .= .-f.P .+ 2 .* f.etac .* (Exxc .+ eta_b .* divV)
            Syy .= .-f.P .+ 2 .* f.etac .* (Eyyc .+ eta_b .* divV)

            Txy .= 2 .* etav   .* Exyv
            Hs  .= 4 .* f.etac .* Eii2

            # ------ Residuals

   @views   f.dVxdtauVx .= (  diff(Txy[2:end-1,:],dims=2)./dy 
                           .+ diff(Sxx,dims=1)./dx)
   @views   f.dVydtauVy .= (  diff(Txy[:,2:end-1],dims=1)./dx 
                           .+ diff(Syy,dims=2)./dy)

            dPdtauP    .= - divV
            dTdtauT    .= (To.-f.T)./dtT .- (diff(f.qx,dims=1)./dx 
                                         .+  diff(f.qy,dims=2)./dy) .+ Hs
            # ------ Updates

            # update with damping
            f.Vx[2:end-1,2:end-1] .+= dtauVx .* (f.dVxdtauVx 
                                      .+ dampx.*f.dVxdtauVx0) 
            f.Vy[2:end-1,2:end-1] .+= dtauVy .* (f.dVydtauVy 
                                      .+ dampy.*f.dVydtauVy0)
            f.P .+= dtauP .* dPdtauP
            f.T .+= dtauT .* dTdtauT

            if (iter % nout ==0) # Check

                err_fu = sqrt(sum(f.dVxdtauVx.^2) + 
                              sum(f.dVydtauVy.^2)) / (nx*(ny-1)+ny*(nx-1))
                err_fp = norm(dPdtauP)/length(dPdtauP)
                err_fT = norm(dTdtauT)/length(dTdtauT)
                err    = [err_fu, err_fp, err_fT]
                push!(errs, [time,err_fu])
                if max(err...) < epsi
                    break
                end

                println("-------------------------")
                @printf(" iter  = %d    \n", iter  )
                @printf(" f_{u} = %1.3e \n", err_fu)
                @printf(" f_{p} = %1.3e \n", err_fp)
                @printf(" f_{T} = %1.3e \n", err_fT)

            end

        end


    end

    errs

end
