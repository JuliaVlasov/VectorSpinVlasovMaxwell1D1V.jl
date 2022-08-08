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

    k = fftfreq(M, M) .* 2π ./ L

    #####################################################
    # use FFT to compute A3_x; A3_xx
    partialA3 = zeros(M)
    partial2A3 = zeros(M)
    partialA3 .= real(ifft(1im .* k .* A3))
    partial2A3 .= real(ifft( - k.^ 2 .* A3))
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
            (0.5 * f0[:, i] .+ 0.5 * sqrt(3) .* f2[:, i]),
            N,
            translatorv1[:, i],
            H,
        )
        u2[:, i] .= translation(
            (0.5 * f0[:, i] .- 0.5 * sqrt(3) .* f2[:, i]),
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

    f0 .= f0t
    f1 .= f1t
    f2 .= f2t
    f3 .= f3t
end

export H2fhOperator

struct H2fhOperator

    adv::Translator
    partial::Vector{ComplexF64}
    v1::Vector{Float64}
    v2::Vector{Float64}
    u1::Matrix{Float64}
    u2::Matrix{Float64}
    f2::Matrix{ComplexF64}

    function H2fhOperator(adv)

        partial = zeros(ComplexF64, adv.mesh.M)
        v1 = zeros(adv.mesh.M)
        v2 = zeros(adv.mesh.M)
        u1 = zeros(adv.mesh.N, adv.mesh.M)
        u2 = zeros(adv.mesh.N, adv.mesh.M)
        f2 = zeros(ComplexF64, adv.mesh.M, adv.mesh.N)

        new(adv, partial, v1, v2, u1, u2, f2)

    end

end


export step!

"""
compute the subsystem H2

$(SIGNATURES)

"""
function step!(f0, f1, f2, f3, E3, A3, op::H2fhOperator, t, h_int)

    k, dv = op.adv.mesh.k, op.adv.mesh.dv

    op.partial .= -k .^ 2 .* A3
    ifft!(op.partial)

    op.v1 .= t .* h_int .* real(op.partial) ./ sqrt(3)
    op.v2 .= -op.v1

    op.u1 .= 0.5 * f0 .+ 0.5 * sqrt(3) .* f2
    op.u2 .= 0.5 * f0 .- 0.5 * sqrt(3) .* f2

    translation!(op.u1, op.adv, op.v1)
    translation!(op.u2, op.adv, op.v2)

    transpose!(op.f2, f2)

    op.partial .= 1im .* k .* A3
    ifft!(op.partial)

    f0 .= op.u1 .+ op.u2
    f2 .= op.u1 ./ sqrt(3) .- op.u2 ./ sqrt(3)
    f1 .= cos.(t .* real(op.partial')) .* f1 .+ sin.(t .* real(op.partial')) .* f3
    f3 .= -sin.(t .* real(op.partial')) .* f1 .+ cos.(t .* real(op.partial')) .* f3

    fft!(op.f2, 1)
    E3 .-= t .* h_int .* (1im .* k) .* vec(sum(op.f2, dims = 2)) .* dv

end
