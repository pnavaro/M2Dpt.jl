include("../src/mesh.jl")
include("../src/fields.jl")
include("physics.jl")
include("numerics.jl")

function loops( nit :: Int64, m :: Mesh, f :: Fields )

    dx     = mesh.dx :: Float64
    dy     = mesh.dy :: Float64

    Vx_exp  = zeros(Float64,(nx+1,ny+2))
    Vy_exp  = zeros(Float64,(nx+2,ny+1))
    etav    = zeros(Float64,(nx+1,ny+1))
    Exyv    = zeros(Float64,(nx+1,ny+1))
    Exyc    = zeros(Float64,(nx,ny))    
    divV    = zeros(Float64,(nx,ny))    
    Exxc    = zeros(Float64,(nx,ny))    
    Eyyc    = zeros(Float64,(nx,ny))    
    Eii2    = zeros(Float64,(nx,ny))    

    for it = 1:nit 

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

            dVxdx = (f.Vx[i+1,j]-f.Vx[i,j])/dx :: Float64
            dVydy = (f.Vy[i,j+1]-f.Vy[i,j])/dy :: Float64

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
            Eii2[i,j] = 0.5*(Exxc[i,j]^2 + Eyyc[i,j]^2) + Exyc[i,j]^2 
        end 

    end

end

function vectorized( nit :: Int64, m :: Mesh, f :: Fields )

    dx     = mesh.dx :: Float64
    dy     = mesh.dy :: Float64

    Vx_exp  = zeros(Float64,(nx+1,ny+2))
    Vy_exp  = zeros(Float64,(nx+2,ny+1))
    etav    = zeros(Float64,(nx+1,ny+1)) 
    Exyv    = zeros(Float64,(nx+1,ny+1)) 
    Exyc    = zeros(Float64,(nx,ny)) 
    divV    = zeros(Float64,(nx,ny)) 
    Exxc    = zeros(Float64,(nx,ny)) 
    Eyyc    = zeros(Float64,(nx,ny)) 
    Eii2    = zeros(Float64,(nx,ny)) 

    for it = 1:nit 

        @views Vx_exp[:,1]       .= f.Vx[:,1]
        @views Vx_exp[:,2:end-1] .= f.Vx
        @views Vx_exp[:,end]     .= f.Vx[:,end]
        @views Vy_exp[1,:]       .= f.Vy[1,:]
        @views Vy_exp[2:end-1,:] .= f.Vy
        @views Vy_exp[end,:]     .= f.Vy[end,:]

               divV .=  diff(f.Vx,dims=1)/dx
               divV .+= diff(f.Vy,dims=2)/dy

               Exxc .= diff(f.Vx,dims=1)/dx
               Exxc .-= 1/2*divV
               Eyyc .= diff(f.Vy,dims=2)/dy
               Eyyc .-= 1/2*divV

               Exyv  .= diff(Vx_exp,dims=2)
               Exyv ./= dy
               Exyv .+= diff(Vy_exp,dims=1)
               Exyv ./= dx
               Exyv .*= 0.5

        @views Exyc  .= 0.25 * ( Exyv[1:end-1,1:end-1] 
                               + Exyv[2:end  ,1:end-1]
                               + Exyv[1:end-1,2:end  ]
                               + Exyv[2:end  ,2:end  ])

               Eii2  .= Exxc.^2
               Eii2 .+= Eyyc.^2
               Eii2 .*= 0.5
               Eii2 .+= Exyc.^2

    end

end

using BenchmarkTools

mesh = Mesh( Lx, nx, Ly, ny)

fields = Fields( mesh, Vbc, r, Tamp )

using InteractiveUtils

@code_warntype vectorized(1, mesh, fields)

@btime vectorized( 1000, mesh, fields )
@btime loops( 1000, mesh, fields, )

