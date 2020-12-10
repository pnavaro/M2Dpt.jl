## This is a trick to make a lazy `diff` which does not allocate in broedcasts:
struct Diff{T,N,DIMS} <: AbstractArray{T,N}
    a::Array{T,N}
end
# Note that using the keyword constructor can lead to allocations in some instances as the
# compiler can then fail to infer DIMS
@inline Diff(a::A; dims=1) where A<:AbstractArray{T,N} where {T,N} = Diff{T,N,dims}(a)
@inline Diff(a::A, dims=1) where A<:AbstractArray{T,N} where {T,N} = Diff{T,N,dims}(a)

# define array interface
# https://docs.julialang.org/en/v1/manual/interfaces/
# follows how views are implemented in Base. E.g.
# https://github.com/JuliaLang/julia/blob/aa2a35a3cadde1f77b8d5816ca0e59c5637c16cc/base/subarray.jl#L263
@inline Base.IndexStyle(::Type{<:Diff}) = IndexCartesian()

@inline function Base.size(D::Diff{T,N,DIMS}) where {T,N,DIMS}
    sz = size(D.a)
    return tupledec(sz, DIMS)
end

@inline function Base.getindex(D::Diff{T,N,DIMS}, I::Vararg{Int, N}) where {T,N,DIMS}
    @boundscheck checkbounds(D, I...)
    a = D.a
    dI = tupleinc(I, DIMS)
    @inbounds r = a[dI...] - a[I...]
    return r
end

"Increment number at index i of tuple"
tupleinc(t::Tuple, i) = tupleadd(t, i, 1)
"Decrement number at index i of tuple"
tupledec(t::Tuple, i) = tupleadd(t, i, -1)

"Add number to tuple at one index"
function tupleadd end
tupleadd(t::NTuple{1}, i, n) = (t[1]+n,)
tupleadd(t::NTuple{2}, i, n) = i==1 ? (t[1]+n, t[2]) : (t[1], t[2]+n)
function tupleadd(t::NTuple{3}, i, n)
    if i==1
        (t[1]+n, t[2], t[3])
    elseif i==2
        (t[1], t[2]+n, t[3])
    else
        (t[1], t[2], t[3]+n)
    end
end
function tupleadd(t::NTuple{4}, i, n)
    if i==1
        (t[1]+n, t[2], t[3], t[4])
    elseif i==2
        (t[1], t[2]+n, t[3], t[4])
    elseif i==3
        (t[1], t[2], t[3]+n, t[4])
    else
        (t[1], t[2], t[3], t[4]+n)
    end
end

@views function vectorized_solve_v2( m :: Mesh, f :: Fields )

    time   = 0.0
    dx, dy = m.dx, m.dy
    nx, ny = m.nx, m.ny

    dtT = tetT*1/4.1*min(dx,dy)^2 # explicit timestep for 2D diffusion

    errs  = []

    etav    = zeros(Float64,(nx+1,ny+1))
    Exyv    = zeros(Float64,(nx+1,ny+1))
    Exyc    = zeros(Float64,(nx,ny))
    divV    = zeros(Float64,(nx,ny))
    Exxc    = zeros(Float64,(nx,ny))
    Eyyc    = zeros(Float64,(nx,ny))
    Eii2    = zeros(Float64,(nx,ny))
    etac_phys = zeros(Float64,(nx,ny))
    Sxx     = zeros(Float64,(nx,ny))
    Syy     = zeros(Float64,(nx,ny))
    Txy     = similar(Exyv)
    Hs      = zeros(Float64,(nx,ny))
    To      = similar(f.T)
    dPdtauP = similar(divV)
    dTdtauT = similar(Hs)
    dtauP   = similar(f.etac)
    dtauT   = similar(f.etac)
    dtauVx  = zeros(Float64,(nx-1,ny))
    dtauVy  = zeros(Float64,(nx,ny-1))

    for it = 1:nt # Physical timesteps

        # temperature from previous step (for backward-Euler integration)
        To .= f.T

        time += dtT # update physical time

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

            divV .= (Diff(f.Vx,1)[:,2:end-1] ./ dx
                  .+ Diff(f.Vy,2)[2:end-1,:] ./ dy)

            Exxc .= Diff(f.Vx,1)[:,2:end-1] ./ dx .- 1/2 .* divV
            Eyyc .= Diff(f.Vy,2)[2:end-1,:] ./ dy .- 1/2 .* divV

            Exyv .= 0.5 .* (Diff(f.Vx,2)./dy .+ Diff(f.Vy,1)./dx)

            Exyc .= 0.25 .* ( Exyv[1:end-1,1:end-1] .+
                              Exyv[2:end  ,1:end-1] .+
                              Exyv[1:end-1,2:end  ] .+
                              Exyv[2:end  ,2:end  ])

            # strain rate invariant
            Eii2 .= 0.5 .* (Exxc.^2 .+ Eyyc.^2) .+ Exyc.^2

            # ------ Rheology
            # physical viscosity
            etac_phys .= Eii2.^mpow .* exp.( -f.T.*(1 ./ (1 .+ f.T./T0)) )

            # numerical shear viscosity
            f.etac .= exp.(rel .* log.(etac_phys) .+ (1-rel) .* log.(f.etac))

            # expand viscosity fom cell centroids to vertices
            fill!(etav,0.0)

            etav[2:end-1,2:end-1] .= 0.25*(f.etac[1:end-1,1:end-1]
                                         + f.etac[2:end,2:end]
                                         + f.etac[1:end-1,2:end]
                                         + f.etac[2:end,1:end-1])

            # this cannot be handled as views
            # etav[:      ,[1, end]] .= etav[:        ,[2, end-1]]
            # etav[[1, end],:      ] .= etav[[2, end-1],:        ]
            # thus break it up (but it does not really impact performance):
            etav[:, 1]   .= etav[:, 2]
            etav[:, end] .= etav[:, end-1]
            etav[1, :]   .= etav[2, :]
            etav[end, :] .= etav[end-1, :]

            # ------ Pseudo-Time steps ------

            dtauP   .= tetp .*  4.1 ./ min(nx,ny) .* f.etac .* (1.0+eta_b)

            dtauVx  .= tetv ./ 4.1 .* (min(dx,dy)^2 ./ (
                       0.5 .*(   f.etac[2:end,:]
                              .+ f.etac[1:end-1,:]) ))./(1+eta_b)

            dtauVy  .= tetv ./ 4.1 .* (min(dx,dy)^2 ./ (
                       0.5 .*(    f.etac[:,2:end]
                              .+  f.etac[:,1:end-1]) ))./(1+eta_b)

            dtauT   .= tetT ./ 4.1 .* min(dx,dy)^2

            # ------ Fluxes

            f.qx[2:end-1,:] .= .- Diff(f.T,1) ./ dx
            f.qy[:,2:end-1] .= .- Diff(f.T,2) ./ dy

            Sxx .= .-f.P .+ 2 .* f.etac .* (Exxc .+ eta_b .* divV)
            Syy .= .-f.P .+ 2 .* f.etac .* (Eyyc .+ eta_b .* divV)

            Txy .= 2 .* etav   .* Exyv
            Hs  .= 4 .* f.etac .* Eii2

            # ------ Residuals

            f.dVxdtauVx .= (  Diff(Txy,2)[2:end-1,:]./dy
                           .+ Diff(Sxx,1)./dx)
            f.dVydtauVy .= (  Diff(Txy,1)[:,2:end-1]./dx
                           .+ Diff(Syy,2)./dy)

            dPdtauP    .= - divV
            dTdtauT    .= (To.-f.T)./dtT .- (Diff(f.qx,1)./dx
                                         .+  Diff(f.qy,2)./dy) .+ Hs
            # ------ Updates

            # update with damping
            # f.Vx[2:end-1,2:end-1] .+= dtauVx .* (f.dVxdtauVx
            #                           .+ dampx.*f.dVxdtauVx0)
            # f.Vy[2:end-1,2:end-1] .+= dtauVy .* (f.dVydtauVy
            #                           .+ dampy.*f.dVydtauVy0)
            f.Vx[2:end-1,2:end-1] .= f.Vx[2:end-1,2:end-1] .+ dtauVx .* (f.dVxdtauVx
                                      .+ dampx.*f.dVxdtauVx0)
            f.Vy[2:end-1,2:end-1] .= f.Vy[2:end-1,2:end-1] .+ dtauVy .* (f.dVydtauVy
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

                # println("-------------------------")
                # @printf(" iter  = %d    \n", iter  )
                # @printf(" f_{u} = %1.3e \n", err_fu)
                # @printf(" f_{p} = %1.3e \n", err_fp)
                # @printf(" f_{T} = %1.3e \n", err_fT)

            end

        end


    end

    errs

end
