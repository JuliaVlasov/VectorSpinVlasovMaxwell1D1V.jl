"""
$(SIGNATURES)
"""
function HAA!(f0, f1, f2, f3, E2, E3, A2, A3, t, L, H)

    N, M = size(f0)

    k = fftfreq(M, M) .* 2pi ./ L
    partialA2 = 1im .* k .* A2
    ifft!(partialA2)
    partialA3 = 1im .* k .* A3
    ifft!(partialA3)
    AA2 = real(ifft(A2))
    AA3 = real(ifft(A3))

    v = real(partialA2) .* AA2 .+ real(partialA3) .* AA3
    v .*= t

    for i = 2:M
        E2[i] += t * k[i]^2 * A2[i]
        E3[i] += t * k[i]^2 * A3[i]
    end

    for i = 1:M
        s = sum(view(f0, :, i))
        AA2[i] = 2H / N * AA2[i] * s
        AA3[i] = 2H / N * AA3[i] * s
    end
    E2 .+= t * fft(AA2)
    E3 .+= t * fft(AA3)

    translation!(f0, v, H)
    translation!(f1, v, H)
    translation!(f2, v, H)
    translation!(f3, v, H)

end


export HAAOperator

struct HAAOperator

    adv::Translator
    A2::Vector{Float64}
    A3::Vector{Float64}
    dA2::Vector{ComplexF64}
    dA3::Vector{ComplexF64}
    delta::Vector{Float64}

    function HAAOperator(adv)

        A2 = zeros(adv.mesh.M)
        A3 = zeros(adv.mesh.M)
        dA2 = zeros(ComplexF64, adv.mesh.M)
        dA3 = zeros(ComplexF64, adv.mesh.M)
        delta = zeros(adv.mesh.M)

        new(adv, A2, A3, dA2, dA3, delta)

    end

end

"""
$(SIGNATURES)
"""
function step!(f0, f1, f2, f3, E2, E3, A2, A3, op::HAAOperator, t)

    M = op.adv.mesh.M
    k = op.adv.mesh.k
    dv = op.adv.mesh.dv

    op.dA2 .= 1im .* k .* A2
    ifft!(op.dA2)
    op.dA3 .= 1im .* k .* A3
    ifft!(op.dA3)
    op.A2 .= real(ifft(A2))
    op.A3 .= real(ifft(A3))

    op.delta .= real(op.dA2) .* op.A2 .+ real(op.dA3) .* op.A3
    op.delta .*= t

    for i = 2:M
        E2[i] += t * k[i]^2 * A2[i]
        E3[i] += t * k[i]^2 * A3[i]
    end

    for i = 1:M
        s = sum(view(f0, :, i))
        op.A2[i] = dv * op.A2[i] * s
        op.A3[i] = dv * op.A3[i] * s
    end

    E2 .+= t * fft(op.A2)
    E3 .+= t * fft(op.A3)

    translation!(f0, op.adv, op.delta)
    translation!(f1, op.adv, op.delta)
    translation!(f2, op.adv, op.delta)
    translation!(f3, op.adv, op.delta)

end
