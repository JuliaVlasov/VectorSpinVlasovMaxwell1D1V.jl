using LinearAlgebra

"""
$(SIGNATURES)

interpolate df(x - delta)
"""
function translation!(df, delta, H)

    # first recover oldvector & get the coefficients of piecewise polynomials

    N, M = size(df)
    @assert length(delta) == M

    a = zeros(N)
    b = zeros(N)
    c = zeros(N + 1)
    cc = zeros(N)
    diagonal = ones(N + 1)
    upright = 1 / 3 .* ones(N)
    diagonal[1] = 2 / 3
    diagonal[N+1] = 2 / 3
    diagonal[2:N] .= 4 / 3
    A = SymTridiagonal(diagonal, upright)

    @inbounds for j = 1:M

        c[2:N] .= view(df, 1:N-1, j) .+ view(df, 2:N, j)
        c[1] = df[1, j]
        c[N+1] = df[N, j]
        c .= A \ c

        for i = 2:N
            cc[i] = (-1)^i * (c[i] - c[i-1])
        end
        for i = 2:N
            b[i] = (-1)^i * 2 * sum(view(cc, 2:i))
        end
        b[1] = 0
        a[1:N-1] .= 1 / 2 .* (view(b, 2:N) .- view(b, 1:N-1))
        a[N] = -1 / 2 * b[N]

        for i = 1:N

            beta = i + delta[j] / (2H / N)
            newbeta = beta - N * floor(Int, beta / N)

            if (abs(newbeta) < 1e-20) || (abs(newbeta - N) < 1e-20)
                df[i, j] = df[N, j]
            elseif newbeta >= 1.0
                l = floor(Int, newbeta)
                k = 1 - (newbeta - l)
                val = a[l] / 3 + b[l] / 2 + c[l]
                val += -a[l] / 3 * (1 - k)^3 - b[l] / 2 * (1 - k)^2
                val += -c[l] * (1 - k)
                val += a[l+1] / 3 * (1 - k)^3 + b[l+1] / 2 * (1 - k)^2
                val += c[l+1] * (1 - k)
                df[i, j] = val

            else
                l = N
                k = 1 - newbeta
                val = a[l] / 3 + b[l] / 2 + c[l]
                val += -a[l] / 3 * (1 - k)^3 - b[l] / 2 * (1 - k)^2
                val += -c[l] * (1 - k)
                val += a[1] / 3 * (1 - k)^3 + b[1] / 2 * (1 - k)^2 + c[1] * (1 - k)
                df[i, j] = val
            end
        end
    end

end

export Translator

struct Translator

    mesh::Mesh
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    diag::Vector{Float64}
    dsup::Vector{Float64}
    A::SymTridiagonal{Float64,Vector{Float64}}

    function Translator(mesh::Mesh)

        a = zeros(mesh.N)
        b = zeros(mesh.N)
        c = zeros(mesh.N + 1)
        d = zeros(mesh.N)

        diag = fill(4 / 3, mesh.N + 1)
        diag[begin] = 2 / 3
        diag[end] = 2 / 3
        dsup = 1 / 3 .* ones(mesh.N)
        A = SymTridiagonal(diag, dsup)

        new(mesh, a, b, c, d, diag, dsup, A)

    end


end

"""
$(SIGNATURES)

interpolate df(x - delta) with Parabolic Spline Method (PSM) 

We consider a
linear advection problem in ``p`` direction
```math
\\frac{\\partial f}{\\partial t} + a \\frac{\\partial f}{\\partial x} =0.
```
From the conservation of the volume, we have the following identity

```math
f_{j,\\ell}(t)=\\frac{1}{\\Delta p} \\int_{p_{\\ell-1/2}} ^{p_{\\ell+1/2}} f(x_j,p,t)\\mathrm{d}{p} =\frac{1}{\\Delta p} \\int_{p_{\\ell-1/2}-at} ^{p_{\\ell+1/2}-at} f(x_j,p,0)\\mathrm{d}{p}.
```

For simplicity, denote by ``q\\in [1, M]`` the index such that
``p_{\\ell+1/2}-at \\in [p_{q-1/2},p_{q+1/2}]`` i.e.
``p_{\\ell+1/2}-at \\in C_q``, then we have

```math
f_{j,\\ell}(t) =\\frac{1}{\\Delta p} \\int_{p_{q-1/2}-at} ^{p_{q-1/2}} f(x_j,p,0)\\mathrm{d}{p}+f_{j,q}(0)-\\frac{1}{\\Delta p} \\int_{p_{q+1/2}} ^{p_{q+1/2}-at} f(x_j,p,0)\\mathrm{d}{p}.
```

Here we need to reconstruct a polynomial function ``f(x_j,p,0)`` using the
averages ``f_{j,l}(0)`` using the PSM approach. 

"""
function translation!(df, adv::Translator, delta)

    N, M = adv.mesh.N, adv.mesh.M
    dv = adv.mesh.dv

    @inbounds for j = 1:M

        adv.c[1] = df[1, j]
        for i = 2:N
            adv.c[i] = df[i-1, j] + df[i, j]
        end
        adv.c[N+1] = df[N, j]
        # get result
        adv.c .= adv.A \ adv.c

        for i = 2:N
            adv.d[i] = (-1)^i * (adv.c[i] - adv.c[i-1])
        end
        for i = 2:N
            adv.b[i] = (-1)^i * 2 * sum(view(adv.d, 2:i))
        end
        adv.b[1] = 0
        adv.a[1:N-1] .= 1 / 2 .* (view(adv.b, 2:N) .- view(adv.b, 1:N-1))
        adv.a[N] = -1 / 2 * adv.b[N]

        for i = 1:N

            beta = i + delta[j] / dv
            newbeta = beta - N * floor(Int, beta / N)

            if (abs(newbeta) < 1e-20) || (abs(newbeta - N) < 1e-20)
                df[i, j] = df[N, j]
            elseif newbeta >= 1.0
                l = floor(Int, newbeta)
                k = 1 - (newbeta - l)
                val = adv.a[l] / 3 + adv.b[l] / 2 + adv.c[l]
                val += -adv.a[l] / 3 * (1 - k)^3 - adv.b[l] / 2 * (1 - k)^2
                val += -adv.c[l] * (1 - k)
                val += adv.a[l+1] / 3 * (1 - k)^3 + adv.b[l+1] / 2 * (1 - k)^2
                val += adv.c[l+1] * (1 - k)
                df[i, j] = val
            else
                l = N
                k = 1 - newbeta
                val = adv.a[l] / 3 + adv.b[l] / 2 + adv.c[l]
                val += -adv.a[l] / 3 * (1 - k)^3 - adv.b[l] / 2 * (1 - k)^2
                val += -adv.c[l] * (1 - k)
                val +=
                    adv.a[1] / 3 * (1 - k)^3 + adv.b[1] / 2 * (1 - k)^2 + adv.c[1] * (1 - k)
                df[i, j] = val
            end
        end
    end

end
