struct Mesh
    
    nx :: Int64
    Lx :: Float64
    dx :: Float64
    xc :: Vector{Float64}
    xn :: Vector{Float64}
    ny :: Int64
    Ly :: Float64
    dy :: Float64
    yc :: Vector{Float64}
    yn :: Vector{Float64}

    function Mesh( Lx :: Float64, nx :: Int64, Ly :: Float64, ny :: Int64 )

        dx = Lx/nx
        dy = Ly/ny
        xn = collect(0:dx:Lx)
        xc = collect(dx/2:dx:Lx-dx/2)
        yn = collect(0:dy:Ly)
        yc = collect(dy/2:dy:Ly-dy/2;)

        new( nx, Lx, dx, xc, xn, ny, Ly, dy, yc, yn)

    end

end
