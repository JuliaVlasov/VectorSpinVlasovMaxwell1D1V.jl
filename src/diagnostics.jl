export diagnostics

function diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, M, N, L, H, h_int)

    value1 = 1:(M-1)÷2 .+ 1
    value2 = (M-1)÷2 .+ 2:M

    B2value = zeros(ComplexF64, M)
    B3value = zeros(ComplexF64, M)

    B2value[value1] .= -((2pi * 1im / L * (value1 .- 1))) .* A3[value1]
    B2value[value2] .= -((2pi * 1im / L * (value2 .- 1 .- M))) .* A3[value2]
    B3value[value1] .= ((2pi * 1im / L * (value1 .- 1))) .* A2[value1]
    B3value[value2] .= ((2pi * 1im / L * (value2 .- 1 .- M))) .* A2[value2]

    EE1value = zeros(M)
    EE2value = zeros(M)
    EE3value = zeros(M)
    BB2value = zeros(M)
    BB3value = zeros(M)
    AA2value = zeros(M)
    AA3value = zeros(M)
    Tvaluet = zeros(M)
    uvaluebar = zeros(M)

    EE1value .= real(ifft(vec(E1)))
    EE2value .= real(ifft(vec(E2)))
    EE3value .= real(ifft(vec(E3)))
    BB2value .= real(ifft(vec(B2value)))
    BB3value .= real(ifft(vec(B3value)))
    AA2value .= real(ifft(vec(A2)))
    AA3value .= real(ifft(vec(A3)))

    # electric energy related to E1
    Ex_energy = 1 / 2 * sum(EE1value[:] .^ 2) * L / M
    # electric energy
    E_energy =
        1 / 2 * sum(EE1value[:] .^ 2) * L / M +
        1 / 2 * sum(EE3value[:] .^ 2) * L / M +
        1 / 2 * sum(EE2value[:] .^ 2) * L / M
    # magnetic energy
    B_energy = 1 / 2 * sum(BB2value[:] .^ 2) * L / M + 1 / 2 * sum(BB3value[:] .^ 2) * L / M
    energy2 = E_energy + B_energy
    v1node = (1:N) .* 2 .* H ./ N .- H .- H ./ N

    ff0value = zeros(N, M)
    ff2value = zeros(N, M)
    ff3value = zeros(N, M)
    for j = 1:M
        Jvaluebar = 0
        rhovaluebar = 0
        for k = 1:N
            ff0value[k, j] =
                1 / 2 * (f0[k, j] * ((v1node[k]^2 + AA2value[j]^2 + AA3value[j]^2))) * L /
                M *
                2 *
                H / N
            ff2value[k, j] = h_int * f2[k, j] * (-BB2value[j]) * L / M * 2 * H / N
            ff3value[k, j] = -h_int * f3[k, j] * (BB3value[j]) * L / M * 2 * H / N
            Jvaluebar = Jvaluebar + f0[k, j] * v1node[k]
            rhovaluebar = rhovaluebar + f0[k, j]
        end
        uvaluebar[j] = Jvaluebar / rhovaluebar
    end
    # temperature
    for j = 1:M
        for k = 1:N
            Tvaluet[j] = Tvaluet[j] + f0[k, j] * ((v1node[k] - uvaluebar[j])^2) * 2 * H / N
        end
    end
    energy1 = sum(ff0value .+ ff2value .+ ff3value)
    # total energy
    energy = energy1 + energy2
    # spectrum
    Sz = sum(f3) * L / M * (2 * H / N)
    return Ex_energy, E_energy, B_energy, energy, Sz, Tvaluet
end
