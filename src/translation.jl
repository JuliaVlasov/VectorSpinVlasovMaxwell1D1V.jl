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
