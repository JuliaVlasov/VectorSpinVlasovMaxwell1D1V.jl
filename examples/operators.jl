using Plots
using FFTW
using MAT
using ProgressMeter
using VectorSpinVlasovMaxwell1D1V
using TimerOutputs

import VectorSpinVlasovMaxwell1D1V: initialfunction
import VectorSpinVlasovMaxwell1D1V: initialfields
import VectorSpinVlasovMaxwell1D1V: diagnostics
import VectorSpinVlasovMaxwell1D1V: H2fh!
import VectorSpinVlasovMaxwell1D1V: He!
import VectorSpinVlasovMaxwell1D1V: HAA!
import VectorSpinVlasovMaxwell1D1V: H3fh!
import VectorSpinVlasovMaxwell1D1V: H1f!

const to = TimerOutput()

function operators()

    T = 50 # 4000  # final time
    nv = 129   # partition of x
    nx = 129   # partition of v
    vmin, vmax = -2.5, 2.5   # v domain size()
    kkk = 1.2231333040331807  #ke
    xmin, xmax = 0, 4pi / kkk  # x domain size()
    h = 0.04 #time step size()
    nsteps = floor(Int, T / h + 1.1) # time step number
    a = 0.02 # 0.001; perturbation coefficient
    h_int = 0.2 # hbar
    k0 = 2.0 * kkk
    ww = sqrt(1.0 + k0^2.0) # w0
    ata = 0.2

    mesh = Mesh(xmin, xmax, nx, vmin, vmax, nv)
    adv = PSMAdvection(mesh)

    E1, E2, E3, A2, A3 = initialfields( mesh, a, ww, kkk, k0)
    f0, f1, f2, f3 = initialfunction(mesh, a, kkk, ata)

    Ex_energy = Float64[]
    E_energy = Float64[]
    B_energy = Float64[]
    energy = Float64[]
    Sz = Float64[]
    Tvalue = Vector{Float64}[]
    time = Float64[]

    results = Diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)

    H2fh = H2fhOperator(adv)
    He = HeOperator(adv)
    HAA = HAAOperator(adv)
    H3fh = H3fhOperator(adv)
    H1f = H1fOperator(adv)

    @showprogress 1 for i = 1:nsteps # Loop over time

        @timeit to "H2fh" step!(f0, f1, f2, f3, E3, A3, H2fh, h / 2, h_int)
        @timeit to "He" step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h / 2)
        @timeit to "HAA" step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h / 2)
        @timeit to "H3fh" step!(f0, f1, f2, f3, E2, A2, H3fh, h / 2, h_int)
        @timeit to "H1f" step!(f0, f1, f2, f3, E1, H1f, h)
        @timeit to "H3fh" step!(f0, f1, f2, f3, E2, A2, H3fh, h / 2, h_int)
        @timeit to "HAA" step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h / 2)
        @timeit to "He" step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h / 2)
        @timeit to "H2fh" step!(f0, f1, f2, f3, E3, A3, H2fh, h / 2, h_int)
        @timeit to "diagnostics" save!(results, i*h, f0, f2, f3, E1, E2, E3, A2, A3)

    end

    results

end

results = operators()

show(to)
plot(results.time, log.(results.Ex_energy), label="julia v2")

vars = matread(joinpath(@__DIR__,"sVMEata0p2.mat"))

plot!(vec(vars["time"]), log.(vec(vars["Ex_energy"])), label="matlab")
