var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = VectorSpinVlasovMaxwell1D1V","category":"page"},{"location":"api/#VectorSpinVlasovMaxwell1D1V","page":"API","title":"VectorSpinVlasovMaxwell1D1V","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Documentation for VectorSpinVlasovMaxwell1D1V.","category":"page"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [VectorSpinVlasovMaxwell1D1V]","category":"page"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.Mesh","page":"API","title":"VectorSpinVlasovMaxwell1D1V.Mesh","text":"struct Mesh\n\nMesh type to store domain parameters\n\nN::Int64\nNumber of points in v\nM::Int64\nNumber of points in x\nH::Float64\nDomain size v ∈ ]-H,+H[\nL::Float64\nDomain size x ∈ [0,L]\nk::Vector{Float64}\nWave number vector to compute derivative with FFTs\ndx::Float64\nSize step along x\ndv::Float64\nSize step along v\nx::Vector{Float64}\npoints along x direction\nv::Vector{Float64}\npoints along v direction\n\n\n\n\n\n","category":"type"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H1f!-NTuple{8, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H1f!","text":"H1f!(f0, f1, f2, f3, E1, t, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H1f-NTuple{10, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H1f","text":"H1f(f0, f1, f2, f3, E1, t, M, N, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H2fh!-NTuple{10, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H2fh!","text":"compute the subsystem H2\n\nH2fh!(f0, f1, f2, f3, E3, A3, t, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H2fh-NTuple{12, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H2fh","text":"compute the subsystem H2\n\nH2fh(f0, f1, f2, f3, E3, A3, t, M, N, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H3fh!-NTuple{10, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H3fh!","text":"H3fh!(f0, f1, f2, f3, E2, A2, t, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H3fh-NTuple{12, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H3fh","text":"H3fh(f0, f1, f2, f3, E2, A2, t, M, N, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.HAA!-NTuple{11, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.HAA!","text":"HAA!(f0, f1, f2, f3, E2, E3, A2, A3, t, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.HAA-NTuple{13, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.HAA","text":"HAA(f0, f1, f2, f3, E2, E3, A2, A3, t, M, N, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.He!-NTuple{11, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.He!","text":"He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, t, H)\n\n\ncompute the first subsystem He()\n\nM is odd number\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.He-NTuple{13, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.He","text":"He(f0, f1, f2, f3, E1, E2, E3, A2, A3, t, M, N, H)\n\n\ncompute the first subsystem He()\n\nM is odd number\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.numeint-Tuple{Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.numeint","text":"computation of integral average in each cell using newton-cotes formula\n\nnumeint(value, N)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.recover-Tuple{Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.recover","text":"recover(f, N)\n\n\nGiven integral average in each cell;this function could compute coefficients a;b;c of piecewise quadratic polynomial using PSM method\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-NTuple{9, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E2, A2, op, t, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any, HeOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, op, t)\n\n\ncompute the first subsystem He()\n\nM is odd number\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, Any, Any, HAAOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E2, E3, A2, A3, op, t)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, H2fhOperator, Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"compute the subsystem H2\n\nstep!(f0, f1, f2, f3, E3, A3, op, t, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, H1fOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E1, op, t)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.translation!-Tuple{Any, Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.translation!","text":"translation!(df, delta, H)\n\n\ninterpolate df(x - delta)\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.translation!-Tuple{Any, Translator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.translation!","text":"translation!(df, adv, delta)\n\n\ninterpolate df(x - delta)\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.translation-NTuple{4, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.translation","text":"translation(oldvector, N, delta, H)\n\n\noldvector is the integral average value in each cell of function f(x) newvector is the integral average value in each cell of function f(x+delta)\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorSpinVlasovMaxwell1D1V","category":"page"},{"location":"#VectorSpinVlasovMaxwell1D1V","page":"Home","title":"VectorSpinVlasovMaxwell1D1V","text":"","category":"section"},{"location":"#Scalar-spin-laser-plasma-model","page":"Home","title":"Scalar spin laser plasma model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Particle distribution function f(x p mathbf s t), xin 0L pin mathbbR  are scalars mathbf s=(s_1s_2s_3) in mathbbR^3, mathbf E = (E_x mathbf E_perp) = (E_x E_y E_z), mathbf A = (A_x mathbf A_perp) = (0 A_y A_z) and mathbf B =nablatimesmathbf  A = (0- partial_xA_z  partial_xA_y).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The scalar spin Vlasov–Maxwell  system introduced in \\cite{crouseilles2021geometric} is given by ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginequationlabeleqreduced\nleft\nbeginaligned\nfracpartial fpartial t + p fracpartial fpartial x +  E_x - mathfrakh s_2 fracpartial^2 A_zpartial x^2 + mathfrakh s_3 fracpartial^2 A_ypartial x^2  - mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial fpartial p   \n hspace3cm+ s_3 fracpartial A_zpartial x + s_2 fracpartial A_ypartial x -s_1 fracpartial A_ypartial x -s_1 fracpartial A_zpartial x  cdot fracpartial fpartial mathbf s = 0\nfracpartial E_xpartial t = -int_mathbbR^4 p f  mathrmdpmathrmdmathrmmathbf s\nfracpartial E_ypartial t = - fracpartial^2 A_ypartial x^2 + A_y int_mathbbR^4  f  mathrmdpmathrmdmathrmmathbf s + mathfrakhint_mathbbR^4 s_3 fracpartial fpartial xmathrmdpmathrmdmathrmmathbf s\nfracpartial E_zpartial t = - fracpartial^2 A_zpartial x^2 + A_z int_mathbbR^4  f  mathrmdpmathrmdmathrmmathbf s - mathfrakhint_mathbbR^4 s_2 fracpartial fpartial xmathrmdpmathrmdmathrmmathbf s\n fracpartial mathbf A_perppartial t = - mathbf E_perp\nfracpartial E_xpartial x = int_mathbbR^4 f mathrmdpmathrmdmathrmmathbf s - 1  text(Poisson equation)\nendaligned\nright\nendequation","category":"page"},{"location":"example/#Example","page":"Example","title":"Example","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"using Plots\nusing FFTW\nusing VectorSpinVlasovMaxwell1D1V\n\nfunction run()\n\n    T = 50 # 4000  # final time\n    M = 65   # partition of x\n    N = 129   # partition of v\n    H = 5.0 / 2   # v domain size()\n    kkk = 1.2231333040331807  #ke\n    L = 4pi / kkk  # x domain size()\n    h = 0.04 #time step size()\n    NUM = floor(Int, T / h + 1.1) # time step number\n    a = 0.02 # 0.001; perturbation coefficient\n    h_int = 0.2 # hbar\n    k0 = 2.0 * kkk\n    ww = sqrt(1.0 + k0^2.0) # w0\n\n    mesh = Mesh(N, M, H, L)\n    adv  = Translator(mesh)\n    \n    E1, E2, E3, A2, A3 = initialfields( mesh, a, ww, kkk, k0)\n    f0, f1, f2, f3 = initialfunction(mesh, a, kkk, ata)\n\n    Ex_energy = Float64[]\n    E_energy = Float64[]\n    B_energy = Float64[]\n    energy = Float64[]\n    Sz = Float64[]\n    Tvalue = Vector{Float64}[]\n    time = Float64[]\n\n    results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)\n    push!(Ex_energy, results[1])\n    push!(E_energy, results[2])\n    push!(B_energy, results[3])\n    push!(energy, results[4])\n    push!(Sz, results[5])\n    push!(Tvalue, results[6])\n    push!(time, 0.0)\n\n    H2fh = H2fhOperator(adv)\n    He = HeOperator(adv)\n    HAA = HAAOperator(adv)\n    H3fh = H3fhOperator(adv)\n    H1f = H1fOperator(adv)\n\n    for i = 1:NUM # Loop over time\n\n        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)\n        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)\n        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E1, H1f, h)\n        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)\n        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)\n        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)\n        \n        results = diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)\n        push!(Ex_energy, results[1])\n        push!(E_energy, results[2])\n        push!(B_energy, results[3])\n        push!(energy, results[4])\n        push!(Sz, results[5])\n        push!(Tvalue, results[6])\n        push!(time, i*h)\n\n    end\n\n    time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue\n\nend","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"time, Ex_energy, E_energy, B_energy, energy, Sz, Tvalue = run()","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(time, Ex_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(time, E_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(time, B_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(time, energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(time, Sz)","category":"page"}]
}
