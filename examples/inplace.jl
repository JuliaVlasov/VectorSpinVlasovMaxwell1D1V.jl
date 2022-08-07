using Plots
using FFTW
using ProgressMeter
using TimerOutputs
using MAT

import VectorSpinVlasovMaxwell1D1V: initialfunction
import VectorSpinVlasovMaxwell1D1V: numeint
import VectorSpinVlasovMaxwell1D1V: diagnostics
import VectorSpinVlasovMaxwell1D1V: H2fh!
import VectorSpinVlasovMaxwell1D1V: He!
import VectorSpinVlasovMaxwell1D1V: HAA!
import VectorSpinVlasovMaxwell1D1V: H3fh!
import VectorSpinVlasovMaxwell1D1V: H1f!

const to = TimerOutput()

function inplace()

    T = 50 # 4000  # final time
    M = 129   # partition of x
    N = 129   # partition of v
    H = 5.0 / 2   # v domain size()
    kkk = 1.2231333040331807  #ke
    L = 4pi / kkk  # x domain size()
    h = 0.04 #time step size()
    NUM = floor(Int, T / h + 1.1) # time step number
    a = 0.02 # 0.001; perturbation coefficient
    h_int = 0.2 # hbar
    k0 = 2.0 * kkk
    ww = sqrt(1.0 + k0^2.0) # w0

    x = (0:(M-1)) .* L ./ M #mesh in x direction
    v = (1:N) .* 2 .* H ./ N .- H #mesh in v direction

    E0 = 0.123 * ww # Eref
    E1 = fft(a ./ kkk .* sin.(kkk .* x))
    E2 = fft(E0 .* cos.(k0 .* x))
    E3 = fft(E0 .* sin.(k0 .* x))
    A2 = -fft(E0 ./ ww .* sin.(k0 .* x))
    A3 = fft(E0 ./ ww .* cos.(k0 .* x))
    ata = 0.2

    kk = 0.17 # v_th
    f(x, v) = (1 / sqrt(2pi) / kk) * exp(-(v^2 / 2 / kk / kk)) * (1 + a * cos(kkk * x))


    f0 = zeros(N, M)
    f1 = zeros(N, M)
    f2 = zeros(N, M)
    f3 = zeros(N, M)

    for k = 1:M

        s = 0

        for i = 1:N

            v1 = v[k] - 2H / N
            v2 = v[k] - 2H / N * 3 / 4
            v3 = v[k] - 2H / N / 2
            v4 = v[k] - 2H / N / 4
            v5 = v[k] 

            y1 = f(x[i], v1)
            y2 = f(x[i], v2)
            y3 = f(x[i], v3)
            y4 = f(x[i], v4)
            y5 = f(x[i], v5)

            s += (7y1 + 32y2 + 12y3 + 32y4 + 7y5) / 90

            f0[i,k] = s

        end

    end

    f3 .= ata ./ 3.0 .* f0

    # test several properties include electric energy; total energy; spectrum etc. save initial data
    Ex_energy = Float64[]
    E_energy = Float64[]
    B_energy = Float64[]
    energy = Float64[]
    Sz = Float64[]
    Tvalue = Vector{Float64}[]
    time = Float64[]

    results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, M, N, L, H, h_int)
    push!(Ex_energy, results[1])
    push!(E_energy, results[2])
    push!(B_energy, results[3])
    push!(energy, results[4])
    push!(Sz, results[5])
    push!(Tvalue, results[6])
    push!(time, 0.0)

    @showprogress 1 for i = 1:NUM # Loop over time

        @timeit to "H2fh" H2fh!(f0, f1, f2, f3, E3, A3, h / 2, L, H, h_int)
        @timeit to "He" He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, H)
        @timeit to "HAA" HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, L, H)
        @timeit to "H3fh" H3fh!(f0, f1, f2, f3, E2, A2, h / 2, L, H, h_int)
        @timeit to "H1f" H1f!(f0, f1, f2, f3, E1, h, L, H)
        @timeit to "H3fh" H3fh!(f0, f1, f2, f3, E2, A2, h / 2, L, H, h_int)
        @timeit to "HAA" HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h / 2, L, H)
        @timeit to "He" He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h / 2, H)
        @timeit to "H2fh" H2fh!(f0, f1, f2, f3, E3, A3, h / 2, L, H, h_int)

        # save properties each time interation
        results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, M, N, L, H, h_int)
        push!(Ex_energy, results[1])
        push!(E_energy, results[2])
        push!(B_energy, results[3])
        push!(energy, results[4])
        push!(Sz, results[5])
        push!(Tvalue, results[6])
        push!(time, i * h)
    end

    time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue

end


time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue = inplace()

show(to)

plot(time, Ex_energy, label="julia v2")

vars = matread(joinpath(@__DIR__,"sVMEata0p2.mat"))

plot!(vec(vars["time"]), vec(vars["Ex_energy"]), label="matlab")
