"""
$(SIGNATURES)
"""
function HAA(f0, f1, f2, f3, E2, E3, A2, A3, t, M, N, L, H)
    # compute the subsystem HA
    savevaluef0 = copy(f0)
    savevaluef1 = copy(f1)
    savevaluef2 = copy(f2)
    savevaluef3 = copy(f3)
    #####################################################
    #use FFT to compute A2_x; A3_x
    partialA2 = zeros(ComplexF64, M)
    partialA3 = zeros(ComplexF64, M)
    value1 = 1:(M-1)÷2+1
    value2 = (M-1)÷2+2:M
    partialA2[value1] .= (((2pi * 1im / L .* (value1 .- 1))) .* A2[value1])
    partialA2[value2] .= (((2pi * 1im / L .* (value2 .- M .- 1))) .* A2[value2])
    partialA2 .= real(ifft(partialA2))
    partialA3[value1] .= (((2pi * 1im / L .* (value1 .- 1))) .* A3[value1])
    partialA3[value2] .= (((2pi * 1im / L .* (value2 .- M .- 1))) .* A3[value2])
    partialA3 .= real(ifft(partialA3))
    AA = real(ifft(A2))
    AAA = real(ifft(A3))
    # solve transport problem in v direction by Semi-Lagrangain method
    translatorv1 = zeros(N, M)
    nonlinear = real(partialA2) .* (real(ifft(A2))) .+ real(partialA3) .* (real(ifft(A3)))
    for i = 1:M
        translatorv1[:, i] .= t * nonlinear[i]
    end

    f0t = zeros(N, M)
    f1t = zeros(N, M)
    f2t = zeros(N, M)
    f3t = zeros(N, M)
    for i = 1:M
        f0t[:, i] .= translation((savevaluef0[:, i]), N, translatorv1[:, i], H)
        f1t[:, i] .= translation((savevaluef1[:, i]), N, translatorv1[:, i], H)
        f2t[:, i] .= translation((savevaluef2[:, i]), N, translatorv1[:, i], H)
        f3t[:, i] .= translation((savevaluef3[:, i]), N, translatorv1[:, i], H)
    end
    #####################################################

    #cpmputation of E2
    E2t = zeros(ComplexF64,M)
    value1 = 2:(M-1)÷2+1
    value2 = (M-1)÷2+2:M
    E2t[value1] .= E2[value1] .+ t * ((2pi / L * (value1 .- 1)) .^ 2) .* A2[value1]
    E2t[value2] .= E2[value2] .+ t * ((2pi / L * (value2 .- M .- 1)) .^ 2) .* A2[value2]
    E2t[1] = E2[1]
    #cpmputation of E3
    E3t = zeros(ComplexF64,M)
    E3t[value1] .= E3[value1] .+ t * ((2pi / L * (value1 .- 1)) .^ 2) .* A3[value1]
    E3t[value2] .= E3[value2] .+ t * ((2pi / L * (value2 .- M .- 1)) .^ 2) .* A3[value2]
    E3t[1] = E3[1]
    ################################

    II = zeros(M)
    for i = 1:M
        II[i] = 2 * H / N * AA[i] * sum(f0[:, i])
    end
    E2t .= E2t .+ t * fft(II)
    ################################
    for i = 1:M
        II[i] = 2 * H / N * AAA[i] * sum(f0[:, i])
    end
    E3t .= E3t .+ t * fft(II)
    return f0t, f1t, f2t, f3t, E2t, E3t
end



"""
$(SIGNATURES)
"""
function HAA!(f0, f1, f2, f3, E2, E3, A2, A3, t, L, H)
   
    N, M = size(f0)

    # compute the subsystem HA
    #####################################################
    #use FFT to compute A2_x; A3_x
    partialA2 = zeros(ComplexF64, M)
    partialA3 = zeros(ComplexF64, M)
    k = fftfreq(M, M) .* 2pi ./ L
    partialA2 .= 1im .* k .* A2
    ifft!(partialA2)
    partialA3 .= 1im .* k .* A3
    ifft!(partialA3)
    AA = real(ifft(A2))
    AAA = real(ifft(A3))
    # solve transport problem in v direction by Semi-Lagrangain method
    nonlinear = real(partialA2) .* (real(ifft(A2))) .+ real(partialA3) .* (real(ifft(A3)))
    v = t .* nonlinear

    #cpmputation of E2
    value1 = 2:(M-1)÷2+1
    value2 = (M-1)÷2+2:M
    E2[value1] .= E2[value1] .+ t * ((2pi / L * (value1 .- 1)) .^ 2) .* A2[value1]
    E2[value2] .= E2[value2] .+ t * ((2pi / L * (value2 .- M .- 1)) .^ 2) .* A2[value2]
    E2[1] = E2[1]
    #cpmputation of E3
    E3[value1] .= E3[value1] .+ t * ((2pi / L * (value1 .- 1)) .^ 2) .* A3[value1]
    E3[value2] .= E3[value2] .+ t * ((2pi / L * (value2 .- M .- 1)) .^ 2) .* A3[value2]
    E3[1] = E3[1]
    ################################

    II = zeros(M)
    for i = 1:M
        II[i] = 2 * H / N * AA[i] * sum(f0[:, i])
    end
    E2 .= E2 .+ t * fft(II)
    ################################
    for i = 1:M
        II[i] = 2 * H / N * AAA[i] * sum(f0[:, i])
    end
    E3 .= E3 .+ t * fft(II)

    translation!(f0, v, H)
    translation!(f1, v, H)
    translation!(f2, v, H)
    translation!(f3, v, H)

end
