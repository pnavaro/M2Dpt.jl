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

function solve_with_loops( m :: Mesh, f :: Fields )

    time   = 0.0
    dx, dy = mesh.dx, mesh.dy

    dtT = tetT*1/4.1*min(dx,dy)^2 # explicit timestep for 2D diffusion
    
    errs  = []

    Vx_exp  = zeros(Float64,(nx+1,ny+2))
    Vy_exp  = zeros(Float64,(nx+2,ny+1))
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

    @show sum(f.T)

    for it = 1:nt # Physical timesteps

        To .= f.T # temperature from previous step (for backward-Euler integration)

        @show time += dtT # update physical time

        for iter = 1:niter # Pseudo-Transient cycles

            # used for damping x momentum residuals
            f.dVxdtauVx0 .= f.dVxdtauVx .+ dampx .* f.dVxdtauVx0  

            # used for damping y momentum residuals
            f.dVydtauVy0 .= f.dVydtauVy .+ dampy .* f.dVydtauVy0 

            #  Kinematics

            for i = 1:nx+1
                Vx_exp[i,1]  = f.Vx[i,1]
                for j = 1:ny
                    Vx_exp[i,j+1] = f.Vx[i,j]
                end
                Vx_exp[i,ny+2] = f.Vx[i,ny]
            end
                
            for j = 1:ny+1
                Vy_exp[1,j] = f.Vy[1,j]
                for i = 1:nx
                    Vy_exp[i+1,j] = f.Vy[i,j]
                end
                Vy_exp[nx+2,j] = f.Vy[nx,j]
            end

            for j = 1:ny, i=1:nx

                dVxdx = (f.Vx[i+1,j]-f.Vx[i,j])/dx
                dVydy = (f.Vy[i,j+1]-f.Vy[i,j])/dy

                divV[i,j] = dVxdx + dVydy

                Exxc[i,j] = dVxdx - 0.5 * (dVxdx + dVydy)

                Eyyc[i,j] = dVydy - 0.5 * (dVxdx + dVydy)
            end

            for j=1:ny+1, i=1:nx+1
                Exyv[i,j] = 0.5*(  (Vx_exp[i,j+1]-Vx_exp[i,j])/dy 
                                 + (Vy_exp[i+1,j]-Vy_exp[i,j])/dx )
            end

            for j=1:ny, i=1:nx
                Exyc[i,j] = 0.25*(Exyv[i,j  ]+Exyv[i+1,j  ]
                                 +Exyv[i,j+1]+Exyv[i+1,j+1])
            end

            for j=1:ny, i=1:nx
                Eii2[i,j] = 0.5*(Exxc[i,j]^2 + Eyyc[i,j]^2) + Exyc[i,j]^2 # strain rate invariant
            end 


            # ------ Rheology
            for j=1:ny, i=1:nx

                 # physical viscosity
                 etac_phys = Eii2[i,j]^mpow*exp( -f.T[i,j]*(1 / (1 + f.T[i,j]/T0)) ) 

                 # numerical shear viscosity
                 f.etac[i,j] = exp(rel*log(etac_phys) + (1-rel)*log(f.etac[i,j]))

            end

            # expand viscosity fom cell centroids to vertices
            fill!(etav,0.0)

            for j = 2:ny, i = 2:nx
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

            dtauP   .= tetp *  4.1 / min(nx,ny)*f.etac*(1.0+eta_b)

            for j = 1:ny, i=1:nx-1
                dtauVx[i,j]  = tetv * 1/4.1 * (min(dx,dy)^2 / ( 
                          0.5*( f.etac[i+1,j] + f.etac[i,j] ) ))/(1+eta_b)
            end

            for j = 1:ny-1, i=1:nx
                dtauVy[i,j] = tetv * 1/4.1 * (min(dx,dy)^2 / (
                           0.5*(  f.etac[i,j+1] + f.etac[i,j]) ))/(1+eta_b)
            end

            dtauT    = tetT * 1/4.1 * min(dx,dy)^2

            # ------ Fluxes

            for j=1:ny, i=2:nx
                f.qx[i,j] = -(f.T[i,j]-f.T[i-1,j])/dx
            end

            for j=2:ny, i=1:nx
                f.qy[i,j] = -(f.T[i,j]-f.T[i,j-1])/dy
            end

            for j=1:ny, i=1:nx
                Sxx[i,j] = -f.P[i,j] + 2 * f.etac[i,j] * (Exxc[i,j] + eta_b*divV[i,j])
                Syy[i,j] = -f.P[i,j] + 2 * f.etac[i,j] * (Eyyc[i,j] + eta_b*divV[i,j])
            end

            for j=1:ny+1, i=1:nx+1
                Txy[i,j] = 2 * etav[i,j]  * Exyv[i,j]
            end

            for j=1:ny, i=1:nx
                Hs[i,j] = 4 * f.etac[i,j] * Eii2[i,j]
            end

            # ------ Residuals

            for j=1:ny, i=1:nx-1
                f.dVxdtauVx[i,j] = (Txy[i+1,j+1]-Txy[i+1,j])/dy + (Sxx[i+1,j]-Sxx[i,j])/dx
            end

            for j=1:ny-1, i=1:nx
                f.dVydtauVy[i,j] = (Txy[i+1,j+1]-Txy[i,j+1])/dx + (Syy[i,j+1]-Syy[i,j])/dy
            end

             
            for j=1:ny, i=1:nx
                dPdtauP[i,j] = - divV[i,j]
                dTdtauT[i,j] = ((To[i,j]-f.T[i,j])/dtT 
                               - ((f.qx[i+1,j]-f.qx[i,j])/dx + (f.qy[i,j+1]-f.qy[i,j])/dy) 
                               + Hs[i,j])
            end
            # ------ Updates

            # update with damping
            f.Vx[2:end-1,:] .+= dtauVx .* (f.dVxdtauVx .+ dampx.*f.dVxdtauVx0) 
            f.Vy[:,2:end-1] .+= dtauVy .* (f.dVydtauVy .+ dampy.*f.dVydtauVy0)
            f.P             .+= dtauP  .* dPdtauP
            f.T             .+= dtauT  .* dTdtauT

            if (iter % nout ==0) # Check

                fu     = hcat(f.dVxdtauVx, transpose(f.dVydtauVy))
                err_fu = norm(fu)/length(fu) 
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
                return

            end

        end
        return


    end

    errs

end

mesh = Mesh( Lx, nx, Ly, ny)

fields = Fields( mesh, Vbc, r, Tamp )

@time solve_with_loops( mesh, fields )

