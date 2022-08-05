"""
compute the subsystem H2

$(SIGNATURES)

"""
function H2fh(f0, f1, f2, f3, E3, A3, t, M, N, L, H, h_int)
    savevaluef0 = copy(f0)
    savevaluef2 = copy(f2)
    #####################################################
    # use FFT to compute A3_x; A3_xx
    partialA3 = zeros(ComplexF64, M)
    partial2A3 = zeros(ComplexF64, M)
    value1 = 1:(M-1)÷2 .+ 1
    value2 = (M-1)÷2 .+ 2:M
    partialA3[value1] .= (((2pi * 1im / L .* (value1 .- 1))) .* A3[value1])
    partialA3[value2] .= (((2pi * 1im / L .* (value2 .- M .- 1))) .* A3[value2])
    partialA3 .= real(ifft(partialA3))
    partial2A3[value1] .= (-((2pi / L * (value1 .- 1)) .^ 2) .* A3[value1])
    partial2A3[value2] .= (-((2pi / L * (value2 .- M .- 1)) .^ 2) .* A3[value2])
    partial2A3 .= real(ifft(partial2A3))
    # solve transport problem in v direction by Semi-Lagrangian method
    translatorv1 = zeros(N, M)
    translatorv2 = zeros(N, M)
    for i = 1:M
        translatorv1[:, i] .= (t * h_int * partial2A3[i] / sqrt(3)) 
        translatorv2[:, i] .= -translatorv1[:, i]
    end

    f0t = zeros(N, M)
    f1t = zeros(N, M)
    f2t = zeros(N, M)
    f3t = zeros(N, M)
    u1 = zeros(N, M)
    u2 = zeros(N, M)

    for i = 1:M
        u1[:, i] .= translation(
            (0.5 * savevaluef0[:, i] .+ 0.5 * sqrt(3) .* savevaluef2[:, i]),
            N,
            translatorv1[:, i],
            H,
        )
        u2[:, i] .= translation(
            (0.5 * savevaluef0[:, i] .- 0.5 * sqrt(3) .* savevaluef2[:, i]),
            N,
            translatorv2[:, i],
            H,
        )
        f0t[:, i] .= u1[:, i] .+ u2[:, i]
        f2t[:, i] .= u1[:, i] ./ sqrt(3) .- u2[:, i] ./ sqrt(3)
        f1t[:, i] .=
            cos.(t * partialA3[i]) * f1[:, i] .+ sin.(t * partialA3[i]) * f3[:, i]
        f3t[:, i] .=
            -sin.(t * partialA3[i]) * f1[:, i] .+ cos.(t * partialA3[i]) * f3[:, i]
    end
    
    ff2 = complex(f2)
    for i = 1:N
        ff2[i, :] .= fft(ff2[i, :])
    end
    #cpmputation of E3
    E3t = zeros(ComplexF64, M)
    E3t[1] = E3[1]
    for i = 2:(M-1)÷2+1
        E3t[i] =
            E3[i] - t * h_int * (1im * 2pi * (i - 1) / L) * sum(ff2[:, i]) * 2 * H / N
    end
    for i = (M-1)÷2+2:M
        k = i - M - 1
        E3t[i] = E3[i] - t * h_int * (1im * 2 * pi * k / L) * sum(ff2[:, i]) * 2 * H / N
    end
    return f0t, f1t, f2t, f3t, E3t
end

"""
compute the subsystem H2

$(SIGNATURES)

"""
function H2fh!(f0, f1, f2, f3, E3, A3, t, L, H, h_int)

    N, M = size(f0)

    #####################################################
    # use FFT to compute A3_x; A3_xx
    partialA3 = zeros(ComplexF64, M)
    partial2A3 = zeros(ComplexF64, M)
    value1 = 1:(M-1)÷2 .+ 1
    value2 = (M-1)÷2 .+ 2:M
    partialA3[value1] .= (((2pi * 1im / L .* (value1 .- 1))) .* A3[value1])
    partialA3[value2] .= (((2pi * 1im / L .* (value2 .- M .- 1))) .* A3[value2])
    partialA3 .= real(ifft(partialA3))
    partial2A3[value1] .= (-((2pi / L * (value1 .- 1)) .^ 2) .* A3[value1])
    partial2A3[value2] .= (-((2pi / L * (value2 .- M .- 1)) .^ 2) .* A3[value2])
    partial2A3 .= real(ifft(partial2A3))
    # solve transport problem in v direction by Semi-Lagrangian method
    translatorv1 = t .* h_int .* real(partial2A3) ./ sqrt(3)
    translatorv2 = -translatorv1

    u1 = 0.5 * f0 .+ 0.5 * sqrt(3) .* f2
    u2 = 0.5 * f0 .- 0.5 * sqrt(3) .* f2

    translation!( u1, translatorv1, H)
    translation!( u2, translatorv2, H)
    f0 .= u1 .+ u2
    f2 .= u1 ./ sqrt(3) .- u2 ./ sqrt(3)
    f1 .= cos.(t .* real(partialA3')) .* f1 .+ sin.(t .* real(partialA3')) .* f3
    f3 .= -sin.(t .* real(partialA3')) .* f1 .+ cos.(t .* real(partialA3')) .* f3
    
    ff2 = complex(f2)
    for i = 1:N
        ff2[i, :] .= fft(ff2[i, :])
    end
    #cpmputation of E3
    for i = 2:(M-1)÷2+1
        E3[i] = E3[i] - t * h_int * (1im * 2pi * (i - 1) / L) * sum(ff2[:, i]) * 2 * H / N
    end
    for i = (M-1)÷2+2:M
        k = i - M - 1
        E3[i] = E3[i] - t * h_int * (1im * 2 * pi * k / L) * sum(ff2[:, i]) * 2 * H / N
    end
end
