using VectorSpinVlasovMaxwell1D1V
using FFTW
using MAT
using Test

import VectorSpinVlasovMaxwell1D1V: initialfunction, numeint
import VectorSpinVlasovMaxwell1D1V: H2fh, H2fh!
import VectorSpinVlasovMaxwell1D1V: He, He!
import VectorSpinVlasovMaxwell1D1V: HAA, HAA!
import VectorSpinVlasovMaxwell1D1V: H3fh, H3fh!
import VectorSpinVlasovMaxwell1D1V: H1f, H1f!

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
const kk = 0.17 # v_th

df(x, v) = exp(-0.5 * v^2 / kk^2) * (1 + a * cos(kkk * x)) / sqrt(2π) / kk

function initialfunction(x, v)

    v1 = v - 2H / N
    v2 = v - H / 2N
    v3 = v
    v4 = v + H / 2N
    v5 = v + H / N

    y1 = df(x, v1)
    y2 = df(x, v2)
    y3 = df(x, v3)
    y4 = df(x, v4)
    y5 = df(x, v5)

    return 7 / 90 * y1 + 16 / 45 * y2 + 2 / 15 * y3 + 16 / 45 * y4 + 7 / 90 * y5

end

@testset "Initialization" begin

    x = (0:(M-1)) .* L ./ M #mesh in x direction
    v1 = (1:N) .* 2 .* H ./ N .- H #mesh in v direction
    tmp = zeros(M)
    E0 = 0.123 * ww # Eref
    E1 = fft(a ./ kkk .* sin.(kkk .* x))
    E2 = fft(E0 .* cos.(k0 .* x))
    E3 = fft(E0 .* sin.(k0 .* x))
    A2 = -fft(E0 ./ ww .* sin.(k0 .* x))
    A3 = fft(E0 ./ ww .* cos.(k0 .* x))

    mesh = Mesh(N, M, H, L)

    @test E0 ≈ 0.123 * ww # Eref
    @test E1 ≈ fft(a ./ kkk .* sin.(kkk .* mesh.x))
    @test E2 ≈ fft(E0 .* cos.(k0 .* mesh.x))
    @test E3 ≈ fft(E0 .* sin.(k0 .* mesh.x))
    @test A2 ≈ -fft(E0 ./ ww .* sin.(k0 .* mesh.x))
    @test A3 ≈ fft(E0 ./ ww .* cos.(k0 .* mesh.x))

    f0 = zeros(N, M)
    for k = 1:M, i = 1:N
        f0[i, k] = initialfunction(k, x, i, v1, kkk, a)
    end

    f0test = similar(f0)
    for k = 1:M, i = 1:N
        f0test[i, k] = f(mesh.x[k], mesh.v[i])
    end

    @test f0test ≈ f0

    # spin related coefficient
    # 5 nodes in each cell in v direction
    v1node = zeros(5N)
    for i = 1:N
        v1node[5*i-4] = v1[i] - 2H / N
        v1node[5*i-3] = v1[i] - (2H / N) * 3 / 4
        v1node[5*i-2] = v1[i] - (2H / N) / 2
        v1node[5*i-1] = v1[i] - (2H / N) / 4
        v1node[5*i] = v1[i]
    end
    # initialize the solution: the interal at each cell in v direction
    f0node = zeros(5N, M)

    for k = 1:M, i = 1:5N
        f0node[i, k] = initialfunction(k, x, i, v1node, kkk, a)
    end
    f0 = zeros(N, M)
    for k = 1:M
        f0[:, k] .= numeint(f0node[:, k], N)
    end




    #f0, f1, f2, f3, E3 = H2fh(f0, f1, f2, f3, E3, A3, h / 2, M, N, L, H, h_int)
    #f0, f1, f2, f3, A2, A3 = He(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, M, N, H)
    #f0, f1, f2, f3, E2, E3 = HAA(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, M, N, L, H)
    #f0, f1, f2, f3, E2 = H3fh(f0, f1, f2, f3, E2, A2, h / 2, M, N, L, H, h_int)
    #f0, f1, f2, f3, E1 = H1f(f0, f1, f2, f3, E1, h, M, N, L, H)
    #f0, f1, f2, f3, E2 = H3fh(f0, f1, f2, f3, E2, A2, h / 2, M, N, L, H, h_int)
    #f0, f1, f2, f3, E2, E3 = HAA(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, M, N, L, H)
    #f0, f1, f2, f3, A2, A3 = He(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, M, N, H)
    #f0, f1, f2, f3, E3 = H2fh(f0, f1, f2, f3, E3, A3, h / 2, M, N, L, H, h_int)

end


@testset "H1f.jl" begin

    mesh = Mesh(N, M, H, L)
    adv = Translator(mesh)

    E0 = 0.123 * ww # Eref
    E1 = fft(a ./ kkk .* sin.(kkk .* mesh.x))
    E2 = fft(E0 .* cos.(k0 .* mesh.x))
    E3 = fft(E0 .* sin.(k0 .* mesh.x))
    A2 = -fft(E0 ./ ww .* sin.(k0 .* mesh.x))
    A3 = fft(E0 ./ ww .* cos.(k0 .* mesh.x))

    fields = matread("fields0.mat")

    @test E1 ≈ fields["E1"]
    @test E2 ≈ fields["E2"]
    @test E3 ≈ fields["E3"]
    @test A2 ≈ fields["A2"]
    @test A3 ≈ fields["A3"]

    f0 = zeros(N, M)
    f1 = zeros(N, M)
    f2 = zeros(N, M)
    f3 = zeros(N, M)

    for k = 1:M, i = 1:N
        f0[i, k] = initialfunction(mesh.x[k], mesh.v[i])
    end

    fields = matread("fields0.mat")

    f3 .= ata ./ 3.0 .* f0

    H2fh = H2fhOperator(adv)
    He = HeOperator(adv)
    HAA = HAAOperator(adv)
    H3fh = H3fhOperator(adv)
    H1f = H1fOperator(adv)

    step!(f0, f1, f2, f3, E3, A3, H2fh, h / 2, h_int)
    step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h / 2)
    step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h / 2)
    step!(f0, f1, f2, f3, E2, A2, H3fh, h / 2, h_int)
    step!(f0, f1, f2, f3, E1, H1f, h)

    H2fh!(f0, f1, f2, f3, E3, A3, h / 2, L, H, h_int)
    He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, H)
    HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, L, H)
    H3fh!(f0, f1, f2, f3, E2, A2, h / 2, L, H, h_int)
    H1f!(f0, f1, f2, f3, E1, h, L, H)
    H3fh!(f0, f1, f2, f3, E2, A2, h / 2, L, H, h_int)
    HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, L, H)
    He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, H)
    H2fh!(f0, f1, f2, f3, E3, A3, h / 2, L, H, h_int)

    @test true

end
