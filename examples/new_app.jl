using Plots
using FFTW
using ProgressMeter

import VectorSpinVlasovMaxwell1D1V: initialfunction
import VectorSpinVlasovMaxwell1D1V: numeint
import VectorSpinVlasovMaxwell1D1V: diagnostics
import VectorSpinVlasovMaxwell1D1V: H2fh!
import VectorSpinVlasovMaxwell1D1V: He!
import VectorSpinVlasovMaxwell1D1V: HAA!
import VectorSpinVlasovMaxwell1D1V: H3fh!
import VectorSpinVlasovMaxwell1D1V: H1f!

function new_main()

    T = 10 # 4000  # final time
    M = 65   # partition of x
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
    v1 = (1:N) .* 2 .* H ./ N .- H #mesh in v direction
    tmp = zeros(M)
    E0 = 0.123 * ww # Eref
    E1 = fft(a ./ kkk .* sin.(kkk .* x))
    E2 = fft(E0 .* cos.(k0 .* x))
    E3 = fft(E0 .* sin.(k0 .* x))
    A2 = -fft(E0 ./ ww .* sin.(k0 .* x))
    A3 = fft(E0 ./ ww .* cos.(k0 .* x))
    ata = 0.2

    # spin related coefficient
    # 5 nodes in each cell in v direction
    v1node = zeros(5N)
    for i = 1:N
        v1node[5*i-4] = v1[i] - 2 * H / N
        v1node[5*i-3] = v1[i] - (2 * H / N) * 3 / 4
        v1node[5*i-2] = v1[i] - (2 * H / N) * 1 / 2
        v1node[5*i-1] = v1[i] - (2 * H / N) * 1 / 4
        v1node[5*i] = v1[i]
    end
    # initialize the solution: the interal at each cell in v direction
    f0_value_at_node = zeros(5N, M)
    f1_value_at_node = zeros(5N, M)
    f2_value_at_node = zeros(5N, M)
    f3_value_at_node = zeros(5N, M)

    for k = 1:M, i = 1:5N
        f0_value_at_node[i, k] = initialfunction(k, x, i, v1node, kkk, a)
        f1_value_at_node[i, k] = 0.0
        f2_value_at_node[i, k] = 0.0
        f3_value_at_node[i, k] = (ata / 3.0) * f0_value_at_node[i, k]
    end

    f0 = zeros(N, M)
    f1 = zeros(N, M)
    f2 = zeros(N, M)
    f3 = zeros(N, M)

    for k = 1:M
        f0[:, k] .= numeint(f0_value_at_node[:, k], N)
        f1[:, k] .= numeint(f1_value_at_node[:, k], N)
        f2[:, k] .= numeint(f2_value_at_node[:, k], N)
        f3[:, k] .= numeint(f3_value_at_node[:, k], N)
    end

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

        H2fh!(f0, f1, f2, f3, E3, A3, h/2, L, H, h_int)
        He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h/2, H)
        HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h/2, L, H)
        H3fh!(f0, f1, f2, f3, E2, A2, h/2, L, H, h_int)
        H1f!(f0, f1, f2, f3, E1, h, L, H)
        H3fh!(f0, f1, f2, f3, E2, A2, h/2, L, H, h_int)
        HAA!(f0, f1, f2, f3, E2, E3, A2, A3, h/2, L, H)
        He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, h/2, H)
        H2fh!(f0, f1, f2, f3, E3, A3, h/2, L, H, h_int)
        
        # save properties each time interation
        results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, M, N, L, H, h_int)
        push!(Ex_energy, results[1])
        push!(E_energy, results[2])
        push!(B_energy, results[3])
        push!(energy, results[4])
        push!(Sz, results[5])
        push!(Tvalue, results[6])
        push!(time, i*h)
    end

    time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue

end


@time time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue = new_main()

plot(time, Ex_energy)
#plot(time, E_energy)
#plot(time, B_energy)
