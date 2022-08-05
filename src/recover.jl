using LinearAlgebra

"""

$(SIGNATURES)

Given integral average in each cell;this function could compute
coefficients a;b;c of piecewise quadratic polynomial using PSM method
"""
function recover(f, N)

    # a;b;c are all row vectors
    f1 = zeros(N + 1)
    f1[2:N] .= f[1:N-1] .+ f[2:N]
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
        b[i] = (-1)^i * 2 * sum(cc[2:i])
    end
    b[1] = 0
    a[1:N-1] .= 1 / 2 .* (b[2:N] .- b[1:N-1])
    a[N] = -1 / 2 * b[N]

    a, b, c
end
