export BSplineAdvection

"""
$(SIGNATURES)

Return the value at x in [0,1] of the B-spline with integer nodes of degree
p with support starting at j.  Implemented recursively using the 
[De Boor's Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)
```math
B_{i,0}(x) := \\left\\{
\\begin{matrix}
1 & \\mathrm{if}  \\quad t_i ≤ x < t_{i+1} \\\\
0 & \\mathrm{otherwise} 
\\end{matrix}
\\right.
```
```math
B_{i,p}(x) := \\frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x) 
+ \\frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```
"""
function bspline(p::Int, j::Int, x::Float64)
    if p == 0
        j == 0 ? (return 1.0) : (return 0.0)
    else
        w = (x - j) / p
        w1 = (x - j - 1) / p
    end
    return w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x)
end

""" 
$(TYPEDEF)

Advection to be computed on each row

$(TYPEDFIELDS)
"""
struct BSplineAdvection <: AbstractAdvection

    mesh::Mesh
    p::Int
    modes::Vector{Float64}
    eig_bspl::Vector{Float64}
    eigalpha::Vector{ComplexF64}
    ft::Matrix{ComplexF64}

    function BSplineAdvection(mesh; p = 3)

        modes = [2pi * i / mesh.nv for i = 0:mesh.nv-1]
        eig_bspl = zeros(mesh.nv)
        eig_bspl .= bspline(p, -div(p + 1, 2), 0.0)
        for i = 1:div(p + 1, 2)-1
            eig_bspl .+= bspline(p, i - (p + 1) ÷ 2, 0.0) * 2 .* cos.(i * modes)
        end
        eigalpha = zeros(ComplexF64, mesh.nv)
        ft = zeros(ComplexF64, mesh.nv, mesh.nx)
        new(mesh, p, modes, eig_bspl, eigalpha, ft)

    end

end

export advection!

"""
$(SIGNATURES)
"""
function advection!(f, adv::BSplineAdvection, edt)

    p = adv.p
    dv :: Float64 = adv.mesh.dv

    adv.ft .= f
    fft!(adv.ft, 1)

    for j in eachindex(edt)

        alpha = edt[j] / dv
        # compute eigenvalues of cubic splines evaluated at displaced points
        ishift = floor(-alpha)
        beta = -ishift - alpha
        fill!(adv.eigalpha, 0.0im)
        for i = -div(p - 1, 2):div(p + 1, 2)
            adv.eigalpha .+= (
                bspline(p, i - div(p + 1, 2), beta) .*
                exp.((ishift + i) .* 1im .* adv.modes)
            )
        end

        # compute interpolating spline using fft and properties of circulant matrices

        adv.ft[:, j] .*= adv.eigalpha ./ adv.eig_bspl

    end

    ifft!(adv.ft, 1)

    f .= real(adv.ft)

end
