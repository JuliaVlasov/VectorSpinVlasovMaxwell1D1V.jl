# Example


```@example test
using Plots
using FFTW
using VectorSpinVlasovMaxwell1D1V

function run()

    T = 50 # 4000  # final time
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

    mesh = Mesh(N, M, H, L)
    adv  = Translator(mesh)
    
    E1, E2, E3, A2, A3 = initialfields( mesh, a, ww, kkk, k0)
    f0, f1, f2, f3 = initialfunction(mesh, a, kkk, ata)

    Ex_energy = Float64[]
    E_energy = Float64[]
    B_energy = Float64[]
    energy = Float64[]
    Sz = Float64[]
    Tvalue = Vector{Float64}[]
    time = Float64[]

    results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)
    push!(Ex_energy, results[1])
    push!(E_energy, results[2])
    push!(B_energy, results[3])
    push!(energy, results[4])
    push!(Sz, results[5])
    push!(Tvalue, results[6])
    push!(time, 0.0)

    H2fh = H2fhOperator(adv)
    He = HeOperator(adv)
    HAA = HAAOperator(adv)
    H3fh = H3fhOperator(adv)
    H1f = H1fOperator(adv)

    for i = 1:NUM # Loop over time

        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)
        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)
        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)
        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)
        step!(f0, f1, f2, f3, E1, H1f, h)
        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)
        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)
        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)
        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)
        
        results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)
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
```

```@example test
time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue = run()
```

```@example test
plot(time, Ex_energy)
```

```@example test
plot(time, E_energy)
```

```@example test
plot(time, B_energy)
```

```@example test
plot(time, energy)
```

```@example test
plot(time, Sz)
```