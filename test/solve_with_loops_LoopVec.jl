using LoopVectorization
using LinearAlgebra

#==
"""
```math
\\frac{∂ υᵢ }{∂xᵢ} = 0, \\
\\frac{τᵢⱼ }{xⱼ} - \\frac{∂ υᵢ }{∂xᵢ}  = 0, \\
τᵢⱼ ϵ̇ᵢⱼ + \\frac{τᵢⱼ }{xⱼ} - \\frac{∂ υᵢ }{∂xᵢ}  = 0, \\
ϵ̇ᵢⱼ  = \\frac{1}{2} \\big( \\frac{∂ υ_i }{∂x_j} + \\frac{∂ υ_j }{∂x_i} \\big))
= 2^{-n} τ\^{n-1}_{II} \\exp \\big( \\frac{n T}{1 + T/T_0} \\big) τ_{ij}
```
"""
==#

function solve_with_loops_v2( m :: Mesh, f :: Fields )

    time   = 0.0
    dx, dy = m.dx, m.dy
    nx, ny = m.nx, m.ny

#    local  dtT :: Float64
    dtT = tetT*1.0/4.1*min(dx,dy)^2 # explicit timestep for 2D diffusion

    etav    = zeros(Float64,(nx+1,ny+1))
    Exyv    = zeros(Float64,(nx+1,ny+1))
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

    @inbounds for it = 1:nt # Physical timesteps

        To .= f.T # temperature from previous step (for backward-Euler integration)
        time += dtT # update physical time

        for iter = 1:niter # Pseudo-Transient cycles

            # used for damping x momentum residuals
            @avx f.dVxdtauVx0 .= f.dVxdtauVx .+ dampx .* f.dVxdtauVx0

            # used for damping y momentum residuals
            @avx f.dVydtauVy0 .= f.dVydtauVy .+ dampy .* f.dVydtauVy0

            #  Kinematics

            for i = 1:nx+1
                f.Vx[i,   1] = f.Vx[i,2]
            end
            for i = 1:nx+1
                f.Vx[i,ny+2] = f.Vx[i,ny+1]
            end
            for j = 1:ny+1
                f.Vy[1,   j] = f.Vy[2,j]
            end
            for j = 1:ny+1
                f.Vy[nx+2,j] = f.Vy[nx+1,j]
            end

            @avx for j = 1:ny, i=1:nx

                dVxdx = (f.Vx[i+1,j+1]-f.Vx[i,j+1])/dx
                dVydy = (f.Vy[i+1,j+1]-f.Vy[i+1,j])/dy

                divV[i,j] = dVxdx + dVydy

                Exxc[i,j] = dVxdx - 0.5 * (dVxdx + dVydy)

                Eyyc[i,j] = dVydy - 0.5 * (dVxdx + dVydy)

            end

            @avx for j=1:ny+1, i=1:nx+1
                Exyv[i,j] = 0.5*(  (f.Vx[i,j+1]-f.Vx[i,j])/dy
                                 + (f.Vy[i+1,j]-f.Vy[i,j])/dx )
            end

            @avx for j=1:ny, i=1:nx
                Exyc = 0.25*(Exyv[i,j  ]+Exyv[i+1,j  ]
                            +Exyv[i,j+1]+Exyv[i+1,j+1])

                # strain rate invariant
                Eii2[i,j] = 0.5*(Exxc[i,j]^2 + Eyyc[i,j]^2) + Exyc^2
            end


            # ------ Rheology
            @avx for i in eachindex(Eii2)

                 # physical viscosity
                 etac_phys = Eii2[i]^mpow*exp(-f.T[i]*(1/(1+f.T[i]/T0)))
                 # numerical shear viscosity
                 f.etac[i] = exp(rel*log(etac_phys)+(1-rel)*log(f.etac[i]))

            end

            # expand viscosity fom cell centroids to vertices
            fill!(etav,0.0)

            @avx for j = 2:ny, i = 2:nx
                etav[i,j] = 0.25*(f.etac[i-1,j-1]+f.etac[i,j-1]
                                 +f.etac[i-1,j  ]+f.etac[i,j  ])
            end

            for i = 1:nx+1
                etav[i,1   ] = etav[i, 2]
                etav[i,ny+1] = etav[i,ny]
            end
            for j = 1:ny+1
                etav[1,   j] = etav[2, j]
                etav[nx+1,j] = etav[nx,j]
            end

            # ------ Pseudo-Time steps ------

            @avx for i in eachindex(dtauP)
                dtauP[i]  = tetp *  4.1 / min(nx,ny)*f.etac[i]*(1.0+eta_b)
            end

            @avx for j = 1:ny, i=1:nx-1
                dtauVx[i,j]  = tetv * 1/4.1 * (min(dx,dy)^2 / (
                          0.5*( f.etac[i+1,j] + f.etac[i,j] ) ))/(1+eta_b)
            end

            @avx for j = 1:ny-1, i=1:nx
                dtauVy[i,j] = tetv * 1/4.1 * (min(dx,dy)^2 / (
                           0.5*(  f.etac[i,j+1] + f.etac[i,j]) ))/(1+eta_b)
            end

            dtauT = tetT * 1/4.1 * min(dx,dy)^2

            # ------ Fluxes

            @avx for j=1:ny, i=2:nx
                f.qx[i,j] = -(f.T[i,j]-f.T[i-1,j])/dx
            end

            @avx for j=2:ny, i=1:nx
                f.qy[i,j] = -(f.T[i,j]-f.T[i,j-1])/dy
            end

            @avx  for i in eachindex(Sxx)
                Sxx[i] = -f.P[i] + 2 * f.etac[i] * (Exxc[i] + eta_b*divV[i])
            end
            @avx  for i in eachindex(Syy)
                Syy[i] = -f.P[i] + 2 * f.etac[i] * (Eyyc[i] + eta_b*divV[i])
            end

            @avx  for i in eachindex(Txy)
                Txy[i] = 2 * etav[i]  * Exyv[i]
            end

            @avx  for i in eachindex(Hs)
                Hs[i] = 4 * f.etac[i] * Eii2[i]
            end

            # ------ Residuals ----------

            @avx for j=1:ny, i=1:nx-1
                f.dVxdtauVx[i,j] = ((Txy[i+1,j+1]-Txy[i+1,j])/dy
                                 + (Sxx[i+1,j]-Sxx[i,j])/dx)
            end

            @avx for j=1:ny-1, i=1:nx
                f.dVydtauVy[i,j] = ((Txy[i+1,j+1]-Txy[i,j+1])/dx
                                 + (Syy[i,j+1]-Syy[i,j])/dy)
            end

            @avx for j=1:ny, i=1:nx
                dPdtauP[i,j] = - divV[i,j]
                dTdtauT[i,j] = ((To[i,j]-f.T[i,j])/dtT
                               - ((f.qx[i+1,j]-f.qx[i,j])/dx
                                + (f.qy[i,j+1]-f.qy[i,j])/dy)
                               + Hs[i,j])
            end
            # ------ Updates ------------

            # update with damping
            @avx for j = 1:ny, i = 1:nx-1
                 f.Vx[i+1,j+1] += dtauVx[i,j] * (f.dVxdtauVx[i,j]
                                         + dampx*f.dVxdtauVx0[i,j])
            end

            @avx for j = 1:ny-1, i = 1:nx
                f.Vy[i+1,j+1] += dtauVy[i,j] * (f.dVydtauVy[i,j] + dampy*f.dVydtauVy0[i,j])
            end

            # f.P .+= dtauP .* dPdtauP
            # f.T .+= dtauT .* dTdtauT

            @avx f.P .= f.P .+ dtauP .* dPdtauP
            @avx f.T .= f.T .+ dtauT .* dTdtauT


            if (iter % nout ==0) # Check

                err_fu = sqrt(sum(f.dVxdtauVx.^2)
                             +sum(f.dVydtauVy.^2))/(nx*(ny-1)+ny*(nx-1))
                err_fp = norm(dPdtauP)/length(dPdtauP)
                err_fT = norm(dTdtauT)/length(dTdtauT)
                err    = [err_fu, err_fp, err_fT]
                if max(err...) < epsi
                    # println("-------------------------")
                    # @printf(" iter  = %d    \n", iter  )
                    # @printf(" f_{u} = %1.3e \n", err_fu)
                    # @printf(" f_{p} = %1.3e \n", err_fp)
                    # @printf(" f_{T} = %1.3e \n", err_fT)
                    break
                end


            end

        end



    end

end
