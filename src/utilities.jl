"""
    expand!(Vexp, V, dims)
   
expand array using BC's - Free slip

"""
function expand!( Vexp::Array{Float64,2}, V::Array{Float64,2}, dims::Int64 )

    if dims == 1
        Vexp[:,2:end-1] .= V
        Vexp[:,1]       .= @view V[:,2]
        Vexp[:,end]     .= @view V[:,end-1]
    else
        Vexp[2:end-1,:] .= V
        Vexp[1,:]       .= @view V[2,:]
        Vexp[end,:]     .= @view V[end-1,:]
    end

end

