
"""
$(SIGNATURES)
"""
function H3fh(f0, f1, f2, f3, E2, A2, t, M, N, L, H, h_int)
    # compute the subsystem H3
    savevaluef0 = copy(f0)
    savevaluef3 = copy(f3)
    #####################################################
    # use FFT to compute A2_x; A2_xx
    partialA2 = zeros(ComplexF64, M)
    partial2A2 = zeros(ComplexF64, M)
    value1 = 1:(M-1)÷2+1
    value2 = (M-1)÷2+2:M
    partialA2[value1] = (((2pi * 1im / L * (value1 .- 1))) .* A2[value1])
    partialA2[value2] = (((2pi * 1im / L * (value2 .- M .- 1))) .* A2[value2])
    partialA2 = real(ifft(partialA2))
    partial2A2[value1] = (-((2pi / L * (value1 .- 1)) .^ 2) .* A2[value1])
    partial2A2[value2] = (-((2pi / L * (value2 .- M .- 1)) .^ 2) .* A2[value2])
    partial2A2 = real(ifft(partial2A2))
    # solve transport problem in v direction by Semi-Lagrangain method
    translatorv1 = zeros(N, M)
    translatorv2 = zeros(N, M)
    for i = 1:M
        translatorv1[:, i] .= -(t * h_int * real(partial2A2[i]) / sqrt(3)) 
        translatorv2[:, i] .= -translatorv1[:, i]
    end

    #####################################################
    f0t = zeros(N, M)
    f1t = zeros(N, M)
    f2t = zeros(N, M)
    f3t = zeros(N, M)
    u1 = zeros(N, M)
    u2 = zeros(N, M)
    for i = 1:M
        u1[:, i] .= translation(
            (0.5*savevaluef0[:, i] .+ 0.5*sqrt(3)*savevaluef3[:, i]),
            N,
            translatorv1[:, i],
            H,
        )
        u2[:, i] .= translation(
            (0.5*savevaluef0[:, i] .- 0.5*sqrt(3)*savevaluef3[:, i]),
            N,
            translatorv2[:, i],
            H,
        )
        f0t[:, i] .= u1[:, i] .+ u2[:, i]
        f3t[:, i] .= u1[:, i] ./ sqrt(3) .- u2[:, i] ./ sqrt(3)
        f1t[:, i] .=
        cos(t * partialA2[i]) .* f1[:, i] .+ sin(t * real(partialA2[i])) .* f2[:, i]
        f2t[:, i] .=
        -sin(t * partialA2[i]) .* f1[:, i] .+ cos(t * real(partialA2[i])) .* f2[:, i]
    end
    #####################################################
    ff3 = zeros(ComplexF64, size(f3))
    for i = 1:N
        ff3[i, :] .= fft(ff3[i, :])
    end
    #computation of E3
    E2t = zeros(ComplexF64, M)
    E2t[1] = E2[1]
    for i = 2:(M-1)÷2+1
        E2t[i] =
            E2[i] + t * h_int * (1im * 2pi * (i - 1) / L) * sum(ff3[:, i]) * 2 * H / N
    end
    for i = (M-1)÷2+2:M
        k = i - M - 1
        E2t[i] = E2[i] + t * h_int * (1im * 2pi * k / L) * sum(ff3[:, i]) * 2 * H / N
    end
    return f0t, f1t, f2t, f3t, E2t
end

"""
$(SIGNATURES)
"""
function H3fh!(f0::Matrix{Float64}, f1::Matrix{Float64}, f2::Matrix{Float64}, f3::Matrix{Float64}, E2, A2, t, L, H, h_int)
   
    N, M = size(f0)

    #####################################################
    # use FFT to compute A2_x; A2_xx
    partialA2 = zeros(ComplexF64, M)
    partial2A2 = zeros(ComplexF64, M)
    value1 = 1:(M-1)÷2+1
    value2 = (M-1)÷2+2:M
    partialA2[value1] = (((2pi * 1im / L * (value1 .- 1))) .* A2[value1])
    partialA2[value2] = (((2pi * 1im / L * (value2 .- M .- 1))) .* A2[value2])
    ifft!(partialA2)
    partial2A2[value1] = (-((2pi / L * (value1 .- 1)) .^ 2) .* A2[value1])
    partial2A2[value2] = (-((2pi / L * (value2 .- M .- 1)) .^ 2) .* A2[value2])
    ifft!(partial2A2)
    # solve transport problem in v direction by Semi-Lagrangain method
    translatorv1 = -t * h_int * real(partial2A2) ./ sqrt(3)
    translatorv2 = -translatorv1

    u1 = 0.5*f0 .+ 0.5*sqrt(3)*f3
    u2 = 0.5*f0 .- 0.5*sqrt(3)*f3

    translation!( u1, translatorv1, H)
    translation!( u2, translatorv2, H)

    f0 .= u1 .+ u2
    f3 .= u1 ./ sqrt(3) .- u2 ./ sqrt(3)
    f1 .= cos.(t * real(partialA2')) .* f1 .+ sin.(t * real(partialA2')) .* f2
    f2 .= -sin.(t * real(partialA2')) .* f1 .+ cos.(t * real(partialA2')) .* f2
    #####################################################
    ff3 = complex(f3)
    fft!(ff3,2)
    #computation of E3
    for i = 2:(M-1)÷2+1
        E2[i] = E2[i] + t * h_int * (1im * 2pi * (i - 1) / L) * sum(ff3[:, i]) * 2 * H / N
    end
    for i = (M-1)÷2+2:M
        k = i - M - 1
        E2[i] = E2[i] + t * h_int * (1im * 2pi * k / L) * sum(ff3[:, i]) * 2 * H / N
    end
end
