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
    partial2A2 = -k .^ 2 .* A2
    ifft!(partial2A2)

    # solve transport problem in v direction by Semi-Lagrangain method
    v1 = -t * h_int * real(partial2A2) ./ sqrt(3)
    v2 = -v1

    u1 = 0.5 * f0 .+ 0.5 * sqrt(3) * f3
    u2 = 0.5 * f0 .- 0.5 * sqrt(3) * f3

    translation!(u1, v1, H)
    translation!(u2, v2, H)

    f0 .= u1 .+ u2
    f3 .= u1 ./ sqrt(3) .- u2 ./ sqrt(3)
    u1 .= cos.(t * real(partialA2')) .* f1 .+ sin.(t * real(partialA2')) .* f2
    u2 .= -sin.(t * real(partialA2')) .* f1 .+ cos.(t * real(partialA2')) .* f2

    f1 .= u1
    f2 .= u2

    ff3 = complex(f3)
    fft!(ff3, 2)

    for i = 2:M
        E2[i] += t * h_int * 1im * k[i] * sum(view(ff3, :, i)) * 2H / N
    end
end


export H3fhOperator

struct H3fhOperator

    adv::Translator
    partial::Vector{ComplexF64}
    v1::Vector{Float64}
    v2::Vector{Float64}
    u1::Matrix{Float64}
    u2::Matrix{Float64}
    f3::Matrix{ComplexF64}

    function H3fhOperator(adv)

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

    op.partial .= -k .^ 2 .* A2
    ifft!(op.partial)

    op.v1 .= -t .* h_int .* real(op.partial) ./ sqrt(3)
    op.v2 .= -op.v1

    op.u1 .= 0.5 .* f0 .+ 0.5 .* sqrt(3) .* f3
    op.u2 .= 0.5 .* f0 .- 0.5 .* sqrt(3) .* f3

    translation!(op.u1, op.adv, op.v1)
    translation!(op.u2, op.adv, op.v2)

    op.partial .= 1im .* k .* A2
    ifft!(op.partial)

    f0 .= op.u1 .+ op.u2
    f3 .= op.u1 ./ sqrt(3) .- op.u2 ./ sqrt(3)
    op.u1 .= cos.(t .* real(op.partial')) .* f1 .+ sin.(t .* real(op.partial')) .* f2
    op.u2 .= -sin.(t .* real(op.partial')) .* f1 .+ cos.(t .* real(op.partial')) .* f2

    f1 .= op.u1
    f2 .= op.u2

    transpose!(op.f3, f3)
    fft!(op.f3, 1)

    for i = 2:M
        E2[i] += t * h_int * 1im * k[i] * sum(view(op.f3, i, :)) * dv
    end
end
