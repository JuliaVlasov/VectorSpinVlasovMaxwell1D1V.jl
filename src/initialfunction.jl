"""
computation of integral average in each cell using newton-cotes formula
$(SIGNATURES)
"""
function numeint(value, N)
    integralvalue = 7 * value[1:5:5N-4]
    integralvalue .+= 32 * value[2:5:5N-3]
    integralvalue .+= 12 * value[3:5:5N-2]
    integralvalue .+= 32 * value[4:5:5N-1]
    integralvalue .+= 7 * value[5:5:5N]

    return integralvalue ./ 90
end

function f(x, v, frequency, a)

    kk = 0.17 # v_th
    value =
        (1 / sqrt(2 * pi) / kk) *
        exp(-v^2 / 2 / kk / kk) *
        (1 + a * cos(frequency * x))

    return value

end

export initialfunction

"""
Initial Gaussian function & add a perturbation in x direction

$(SIGNATURES)
"""

function initialfunction(H, L, N, M, a, frequency, ata)


    x = (0:(M-1)) .* L ./ M 
    v = (1:N) .* 2 .* H ./ N .- H 

    # spin related coefficient
    # 5 nodes in each cell in v direction
    vnode = zeros(5N)
    for i = 1:N
        vnode[5*i-4] = v[i] - 2H / N
        vnode[5*i-3] = v[i] - (2H / N) * 3 / 4
        vnode[5*i-2] = v[i] - (2H / N) / 2
        vnode[5*i-1] = v[i] - (2H / N) / 4
        vnode[5*i] = v[i]
    end
    # initialize the solution: the interal at each cell in v direction
    f0node = zeros(5N, M)
    f1node = zeros(5N, M)
    f2node = zeros(5N, M)
    f3node = zeros(5N, M)

    for k = 1:M, i = 1:5N
        f0node[i, k] = f(x[k], vnode[i], frequency, a)
        f1node[i, k] = 0.0
        f2node[i, k] = 0.0
        f3node[i, k] = (ata / 3.0) * f0node[i, k]
    end

    f0 = zeros(N, M)
    f1 = zeros(N, M)
    f2 = zeros(N, M)
    f3 = zeros(N, M)

    for k = 1:M
        f0[:, k] .= numeint(f0node[:, k], N)
        f1[:, k] .= numeint(f1node[:, k], N)
        f2[:, k] .= numeint(f2node[:, k], N)
        f3[:, k] .= numeint(f3node[:, k], N)
    end

    return f0, f1, f2, f3

end

function initialfunction(mesh, a, frequency, ata)

    xmin, xmax = mesh.xmin, mesh.xmax
    vmin, vmax = mesh.vmin, mesh.vmax
    nx, nv = mesh.nx, mesh.nv
    dv = mesh.dv

    f0 = zeros(nv, nx)
    f1 = zeros(nv, nx)
    f2 = zeros(nv, nx)
    f3 = zeros(nv, nx)

    for k = 1:nx
        for i = 1:nv
            v1 = mesh.v[i] - dv
            v2 = mesh.v[i] - dv * 0.75
            v3 = mesh.v[i] - dv * 0.50
            v4 = mesh.v[i] - dv * 0.25
            v5 = mesh.v[i] 

            y1 = f(mesh.x[k], v1, frequency, a)
            y2 = f(mesh.x[k], v2, frequency, a)
            y3 = f(mesh.x[k], v3, frequency, a)
            y4 = f(mesh.x[k], v4, frequency, a)
            y5 = f(mesh.x[k], v5, frequency, a)

            f0[i,k] = (7y1 + 32y2 + 12y3 + 32y4 + 7y5) / 90
        end
    end

    f3 .= (ata / 3.0) .* f0

    return f0, f1, f2, f3

end
