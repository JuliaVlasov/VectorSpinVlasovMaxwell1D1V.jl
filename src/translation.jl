using LinearAlgebra

"""

$(SIGNATURES)

Given integral average in each cell;this function could compute
coefficients a;b;c of piecewise quadratic polynomial using PSM method
"""
function recover(f, N)

    # a;b;c are all row vectors
    f1 = zeros(N + 1)
    f1[2:N] .= view(f, 1:N-1) .+ view(f, 2:N)
    f1[1] = f[1]
    f1[N+1] = f[N]
    diagonal = ones(N + 1)
    upright = 1 / 3 .* ones(N)
    diagonal[1] = 2 / 3
    diagonal[N+1] = 2 / 3
    diagonal[2:N] .= 4 / 3
    # here we use MATLAB function sparse to get matrix A
    A = SymTridiagonal(diagonal, upright)
    # get result
    c = A \ f1

    a = zeros(N)
    b = zeros(N)
    cc = zeros(N)
    for i = 2:N
        cc[i] = (-1)^i * (c[i] - c[i-1])
    end
    for i = 2:N
        b[i] = (-1)^i * 2 * sum(view(cc, 2:i))
    end
    b[1] = 0
    a[1:N-1] .= 1 / 2 .* (view(b, 2:N) .- view(b, 1:N-1))
    a[N] = -1 / 2 * b[N]

    a, b, c
end

"""
$(SIGNATURES)

oldvector is the integral average value in each cell of function f(x)
newvector is the integral average value in each cell of function f(x+delta)
"""
function translation(oldvector, N, delta, H)

    # first recover oldvector & get the coefficients of piecewise polynomials
    a, b, c = recover(oldvector, N)

    newvector = zeros(N)

    for i = 1:N
        beta = i + delta[i] / (2 * H / N)
        newbeta = beta - N * floor(Int, beta / N)

        if (abs(newbeta) < 1e-20) || (abs(newbeta - N) < 1e-20)
            newvector[i] = oldvector[N]
        elseif newbeta >= 1.0
            l = floor(Int, newbeta)
            k = 1 - (newbeta - l)
            valueI = a[l] / 3 + b[l] / 2 + c[l]
            valueI = valueI - a[l] / 3 * (1 - k)^3 - b[l] / 2 * (1 - k)^2
            valueI = valueI - c[l] * (1 - k)
            valueII = a[l+1] / 3 * (1 - k)^3 + b[l+1] / 2 * (1 - k)^2
            valueII = valueII + c[l+1] * (1 - k)
            newvector[i] = valueI + valueII

        else
            l = N
            k = 1 - newbeta
            valueI = a[l] / 3 + b[l] / 2 + c[l]
            valueI = valueI - a[l] / 3 * (1 - k)^3 - b[l] / 2 * (1 - k)^2
            valueI = valueI - c[l] * (1 - k)
            valueII = a[1] / 3 * (1 - k)^3 + b[1] / 2 * (1 - k)^2 + c[1] * (1 - k)
            newvector[i] = valueI + valueII
        end

    end

    return newvector

end

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
        # get result
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

interpolate df(x - delta)
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
