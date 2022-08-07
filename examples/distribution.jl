using VectorSpinVlasovMaxwell1D1V
using Plots
using Test

const M = 129   # partition of x
const N = 129   # partition of v
const H = 5.0 / 2   # v domain size()
const kkk = 1.2231333040331807  #ke
const L = 4pi / kkk  # x domain size()
const h = 0.04 #time step size()
const a = 0.02 # 0.001; perturbation coefficient
const h_int = 0.2 # hbar
const k0 = 2.0 * kkk
const ww = sqrt(1.0 + k0^2.0) # w0
const ata = 0.2
const kk = 0.17

function initialfunction(x, v, k, a)


    v1 = v - H / N
    v2 = v - H / 2N
    v3 = v
    v4 = v + H / 2N
    v5 = v + H / N

    y1 = f(x, v1)
    y2 = f(x, v2)
    y3 = f(x, v3)
    y4 = f(x, v4)
    y5 = f(x, v5)

    return 7 / 90 * y1 + 16 / 45 * y2 + 2 / 15 * y3 + 16 / 45 * y4 + 7 / 90 * y5

end


function initialfunction(k, x, i, v1int, frequency, a)

    value =
        (1 / sqrt(2 * pi) / kk) *
        exp(-(v1int[i])^2 / 2 / kk / kk) *
        (1 + a * cos(frequency * x[k]))

    return value

end

function initialfunction(x, v)

    exp(- 0.5 * v^2 / kk^2) * (1 + a * cos(kkk * x)) / sqrt(2π) / kk

end


dx = L / M
dv = 2H / N
x1 = (0:(M-1)) .* dx 
v1 = (1:N) .* dv .- H #mesh in v direction

mesh = Mesh(N, M, H, L)

@test x1 ≈ mesh.x 
@test v1 ≈ mesh.v 

f0 = zeros(N, M)

# spin related coefficient
# 5 nodes in each cell in v direction
v1node = zeros(5N)
for i = 1:N
    v1node[5*i-4] = v1[i] - dv
    v1node[5*i-3] = v1[i] - dv * 3 / 4
    v1node[5*i-2] = v1[i] - dv * 1 / 2
    v1node[5*i-1] = v1[i] - dv * 1 / 4
    v1node[5*i] = v1[i]
end
# initialize the solution: the interal at each cell in v direction
f0 = zeros(N, M)
f0test = zeros(N, M)

for j = 1:M, i = 1:N
    f0[i, j] = initialfunction(j, x1, i, v1, kkk, a)
    f0test[i, j] = f(x1[j], v1[i])
end

@test f0test ≈ f0 
contourf(x1, v1, f0)

