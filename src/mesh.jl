export Mesh

""" 
$(TYPEDEF)
Mesh type to store domain parameters
$(TYPEDFIELDS)
"""
struct Mesh

    "Number of points in v"
    nv::Int64
    "Number of points in x"
    nx::Int64
    "Domain size v ∈ ]vmin,vmax["
    vmin::Float64
    "Domain size v ∈ ]vmin,vmax["
    vmax::Float64
    "Domain size x ∈ [xmin,xmax]"
    xmin::Float64
    "Domain size x ∈ [xmin,xmax]"
    xmax::Float64
    "Wave number vector to compute derivative with FFTs"
    kx::Vector{Float64}
    "Size step along x"
    dx::Float64
    "Size step along v"
    dv::Float64
    "points along x direction"
    x::Vector{Float64}
    "points along v direction"
    v::Vector{Float64}

    function Mesh(xmin, xmax, nx, vmin, vmax, nv)

        dx = (xmax - xmin) / nx
        dv = (vmax - vmin) / nv
        kx = collect(2π ./ (xmax - xmin) .* fftfreq(nx, nx))
        x = collect(0:(nx-1)) .* (xmax - xmin) ./ nx 
        v = collect(1:nv) .* dv .- 0.5 * (vmax - vmin)

        new(nv, nx, vmin, vmax, xmin, xmax, kx, dx, dv, x, v)

    end
end
