using Plots
using FFTW
using ProgressMeter
using TimerOutputs
using MAT

import VectorSpinVlasovMaxwell1D1V: initialfields
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
    ata = 0.2

    E1, E2, E3, A2, A3 = initialfields( H, L, N, M, a, ww, kkk, k0)
    f0, f1, f2, f3 = initialfunction(H, L, N, M, a, kkk, ata)

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

plot(time, Ex_energy, label="julia v1")

vars = matread(joinpath(@__DIR__,"sVMEata0p2.mat"))

plot!(vec(vars["time"]), vec(vars["Ex_energy"]), label="matlab")
