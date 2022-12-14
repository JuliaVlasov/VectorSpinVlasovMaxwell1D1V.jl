"""
compute the subsystem H2

$(SIGNATURES)

"""
function H2fh!(f0, f1, f2, f3, E3, A3, t, L, H, h_int)

    N, M = size(f0)

    k = fftfreq(M, M) .* 2π ./ L

    partialA3 = real(ifft(1im .* k .* A3))
    partial2A3 = real(ifft( - k.^ 2 .* A3))
    
    v1 = t .* h_int .* partial2A3 ./ sqrt(3)
    v2 = - v1

    u1 = 0.5 .* f0 .+ 0.5 * sqrt(3) .* f2
    u2 = 0.5 .* f0 .- 0.5 * sqrt(3) .* f2

    translation!( u1, v1, H)
    translation!( u2, v2, H)

    ff2 = complex(f2)
    fft!( ff2, 2)
    @inbounds for i = 2:M
        E3[i] += - t * h_int * 1im * k[i] * sum(ff2[:, i]) * 2 * H / N
    end

    f0 .= u1 .+ u2
    f2 .= u1 ./ sqrt(3) .- u2 ./ sqrt(3)

    u1 .= cos.(t * partialA3') .* f1 .+ sin.(t * partialA3') .* f3
    u2 .= -sin.(t * partialA3') .* f1 .+ cos.(t * partialA3') .* f3

    f1 .= u1
    f3 .= u2
end

export H2fhOperator

struct H2fhOperator

    adv::AbstractAdvection
    partial::Vector{ComplexF64}
    v1::Vector{Float64}
    v2::Vector{Float64}
    u1::Matrix{Float64}
    u2::Matrix{Float64}
    f2::Matrix{ComplexF64}

    function H2fhOperator(adv)

        partial = zeros(ComplexF64, adv.mesh.nx)
        v1 = zeros(adv.mesh.nx)
        v2 = zeros(adv.mesh.nx)
        u1 = zeros(adv.mesh.nv, adv.mesh.nx)
        u2 = zeros(adv.mesh.nv, adv.mesh.nx)
        f2 = zeros(ComplexF64, adv.mesh.nx, adv.mesh.nv)

        new(adv, partial, v1, v2, u1, u2, f2)

    end

end


export step!

"""
compute the subsystem H2

$(SIGNATURES)

"""
function step!(f0, f1, f2, f3, E3, A3, op::H2fhOperator, dt, h_int)

    kx :: Vector{Float64} = op.adv.mesh.kx
    dv :: Float64 = op.adv.mesh.dv

    op.partial .= -kx .^ 2 .* A3
    ifft!(op.partial)

    op.v1 .= - h_int .* real(op.partial) ./ sqrt(3)
    op.v2 .= -op.v1

    op.u1 .= 0.5 * f0 .+ 0.5 * sqrt(3) .* f2
    op.u2 .= 0.5 * f0 .- 0.5 * sqrt(3) .* f2

    advection!(op.u1, op.adv, op.v1, dt)
    advection!(op.u2, op.adv, op.v2, dt)

    transpose!(op.f2, f2)

    op.partial .= 1im .* kx .* A3
    ifft!(op.partial)

    f0 .= op.u1 .+ op.u2
    f2 .= op.u1 ./ sqrt(3) .- op.u2 ./ sqrt(3)
    op.u1 .= cos.(dt .* real(op.partial')) .* f1 .+ sin.(dt .* real(op.partial')) .* f3
    op.u2 .= -sin.(dt .* real(op.partial')) .* f1 .+ cos.(dt .* real(op.partial')) .* f3

    fft!(op.f2, 1)
    E3 .-= dt .* h_int .* (1im .* kx) .* vec(sum(op.f2, dims = 2)) .* dv

    f1 .= op.u1
    f3 .= op.u2

end
