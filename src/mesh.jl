export Mesh

""" 
$(TYPEDEF)
Mesh type to store domain parameters
$(TYPEDFIELDS)
"""
struct Mesh

    "Number of points in v"
    N :: Int64
    "Number of points in x"
    M :: Int64
    "Domain size v ∈ ]-H,+H["
    H :: Float64
    "Domain size x ∈ [0,L]"
    L :: Float64
    "Wave number vector to compute derivative with FFTs"
    k :: Vector{Float64}
    "Size step along x"
    dx :: Float64
    "Size step along v"
    dv :: Float64
    "points along x direction"
    x :: Vector{Float64}
    "points along v direction"
    v :: Vector{Float64}

    function Mesh( N, M, H, L)

        dx = L / M
        dv = 2H / N
        k = collect(2π ./ L .* fftfreq(M, M))
        x = collect((0:(M-1)) .* dx )
        v = collect((1:N) .* dv .- H )

        new( N, M, H, L, k, dx, dv, x, v )

    end
end




