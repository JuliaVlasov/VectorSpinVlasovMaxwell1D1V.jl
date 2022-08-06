
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
function H3fh!(f0, f1, f2, f3, E2, A2, t, L, H, h_int)
   
    N, M = size(f0)

    #####################################################
    # use FFT to compute A2_x; A2_xx
    k = 2π ./ L .* fftfreq(M, M)
    partialA2 = 1im .* k .* A2
    ifft!(partialA2)
    partial2A2 = - k .^ 2 .* A2
    ifft!(partial2A2)
   
    # solve transport problem in v direction by Semi-Lagrangain method
    v1 = -t * h_int * real(partial2A2) ./ sqrt(3)
    v2 = -v1

    u1 = 0.5*f0 .+ 0.5*sqrt(3)*f3
    u2 = 0.5*f0 .- 0.5*sqrt(3)*f3

    translation!( u1, v1, H)
    translation!( u2, v2, H)

    ff3 = complex(f3)

    f0 .= u1 .+ u2
    f3 .= u1 ./ sqrt(3) .- u2 ./ sqrt(3)
    f1 .= cos.(t * real(partialA2')) .* f1 .+ sin.(t * real(partialA2')) .* f2
    f2 .= -sin.(t * real(partialA2')) .* f1 .+ cos.(t * real(partialA2')) .* f2
    
    fft!(ff3,2)
    
    for i = 2:M
        E2[i] += t * h_int * 1im * k[i] * sum(view(ff3,:, i)) * 2H / N
    end
end


export H3fhOperator

struct H3fhOperator

    adv :: Translator
    partial :: Vector{ComplexF64}
    v1 :: Vector{Float64}
    v2 :: Vector{Float64}
    u1 :: Matrix{Float64}
    u2 :: Matrix{Float64}
    f3 :: Matrix{ComplexF64}

    function H3fhOperator( adv )

        N, M = adv.mesh.N, adv.mesh.M
        partial = zeros(ComplexF64, M)
        v1 = zeros(M)
        v2 = zeros(M)
        u1 = zeros(N, M)
        u2 = zeros(N, M)
        f3 = zeros(ComplexF64, M, N)

        new(adv, partial, v1, v2, u1, u2, f3)

    end


end




"""
$(SIGNATURES)
"""
function step!(f0, f1, f2, f3, E2, A2, op, t, h_int)
   
    M = op.adv.mesh.M
    dv = op.adv.mesh.dv
    k = op.adv.mesh.k

    op.partial .= - k .^ 2 .* A2
    ifft!(op.partial)
   
    op.v1 .= -t .* h_int .* real(op.partial) ./ sqrt(3)
    op.v2 .= -op.v1

    op.u1 .= 0.5 .* f0 .+ 0.5 .* sqrt(3) .* f3
    op.u2 .= 0.5 .* f0 .- 0.5 .* sqrt(3) .* f3

    translation!( op.u1, op.adv, op.v1)
    translation!( op.u2, op.adv, op.v2)

    op.partial .= 1im .* k .* A2
    ifft!(op.partial)

    transpose!(op.f3, f3)

    f0 .= op.u1 .+ op.u2
    f3 .= op.u1 ./ sqrt(3) .- op.u2 ./ sqrt(3)
    f1 .= cos.(t .* real(op.partial')) .* f1 .+ sin.(t .* real(op.partial')) .* f2
    f2 .= -sin.(t .* real(op.partial')) .* f1 .+ cos.(t .* real(op.partial')) .* f2
    
    fft!(op.f3,1)
    
    for i = 2:M
        E2[i] += t * h_int * 1im * k[i] * sum(view(op.f3,i, :)) * dv
    end
end
