"""
$(SIGNATURES)

old is the integral average value in each cell of function f(x)
new is the integral average value in each cell of function f(x+delta)
"""
function translation!(df, delta, H)

    # first recover oldvector & get the coefficients of piecewise polynomials

    N, M = size(df)
    @assert length(delta) == M

    f1 = zeros(N + 1)
    a = zeros(N)
    b = zeros(N)
    c = zeros(N+1)
    cc = zeros(N)
    diagonal = ones(N + 1)
    upright = 1 / 3 .* ones(N)
    diagonal[1] = 2 / 3
    diagonal[N+1] = 2 / 3
    diagonal[2:N] .= 4 / 3 
    A = SymTridiagonal(diagonal, upright)

    @inbounds for j = 1:M

        f1[2:N] .= df[1:N-1,j] .+ df[2:N, j]
        f1[1] = df[1,j]
        f1[N+1] = df[N,j]
        # get result
        c .= A \ f1

        for i = 2:N
            cc[i] = (-1)^i * (c[i] - c[i-1])
        end
        for i = 2:N
            b[i] = (-1)^i * 2 * sum(cc[2:i])
        end
        b[1] = 0
        a[1:N-1] .= 1 / 2 .* (b[2:N] .- b[1:N-1])
        a[N] = -1 / 2 * b[N]

        for i = 1:N

            beta = i + delta[j] / (2 * H / N)
            loopnumber = floor(Int, beta / N)
            newbeta = beta - N * loopnumber

            if (abs(newbeta) < 1e-20) || (abs(newbeta - N) < 1e-20)
                df[i,j] = df[N,j]
            elseif newbeta >= 1.0
                index = floor(Int, newbeta)
                k = 1 - (newbeta - index)
                valueI = a[index] / 3 + b[index] / 2 + c[index]
                valueI = valueI - a[index] / 3 * (1 - k)^3 - b[index] / 2 * (1 - k)^2
                valueI = valueI - c[index] * (1 - k)
                valueII = a[index+1] / 3 * (1 - k)^3 + b[index+1] / 2 * (1 - k)^2
                valueII = valueII + c[index+1] * (1 - k)
                df[i,j] = valueI + valueII

            else
                index = N
                k = 1 - newbeta
                valueI = a[index] / 3 + b[index] / 2 + c[index]
                valueI = valueI - a[index] / 3 * (1 - k)^3 - b[index] / 2 * (1 - k)^2
                valueI = valueI - c[index] * (1 - k)
                valueII = a[1] / 3 * (1 - k)^3 + b[1] / 2 * (1 - k)^2 + c[1] * (1 - k)
                df[i,j] = valueI + valueII
            end
        end
    end

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
        loopnumber = floor(Int, beta / N)
        newbeta = beta - N * loopnumber

        if (abs(newbeta) < 1e-20) || (abs(newbeta - N) < 1e-20)
            newvector[i] = oldvector[N]
        elseif newbeta >= 1.0
            index = floor(Int, newbeta)
            k = 1 - (newbeta - index)
            valueI = a[index] / 3 + b[index] / 2 + c[index]
            valueI = valueI - a[index] / 3 * (1 - k)^3 - b[index] / 2 * (1 - k)^2
            valueI = valueI - c[index] * (1 - k)
            valueII = a[index+1] / 3 * (1 - k)^3 + b[index+1] / 2 * (1 - k)^2
            valueII = valueII + c[index+1] * (1 - k)
            newvector[i] = valueI + valueII

        else
            index = N
            k = 1 - newbeta
            valueI = a[index] / 3 + b[index] / 2 + c[index]
            valueI = valueI - a[index] / 3 * (1 - k)^3 - b[index] / 2 * (1 - k)^2
            valueI = valueI - c[index] * (1 - k)
            valueII = a[1] / 3 * (1 - k)^3 + b[1] / 2 * (1 - k)^2 + c[1] * (1 - k)
            newvector[i] = valueI + valueII
        end

    end

    return newvector

end
