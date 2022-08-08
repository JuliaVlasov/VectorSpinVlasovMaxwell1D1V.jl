"""
$(SIGNATURES)
"""
function H1f(f0, f1, f2, f3, E1, t, M, N, L, H)
    #####################################################
    # computation of fO & f
    ff0 = complex(f0)
    ff1 = complex(f1)
    ff2 = complex(f2)
    ff3 = complex(f3)
    v = (1:N) .* 2 .* H ./ N .- H .- H ./ N
    # Fourier transform
    for i = 1:N
        ff0[i, :] .= fft(ff0[i, :])
        ff1[i, :] .= fft(ff1[i, :])
        ff2[i, :] .= fft(ff2[i, :])
        ff3[i, :] .= fft(ff3[i, :])
    end
    f0t = zeros(ComplexF64, N, M)
    f1t = zeros(ComplexF64, N, M)
    f2t = zeros(ComplexF64, N, M)
    f3t = zeros(ComplexF64, N, M)
    # frequency: k_fre=[0:M/2-1,-M/2:-1]
    # M is odd; odd frequency number
    k_fre = [0:(M-1)÷2; -(M - 1)÷2:-1]
    # FFT in x direction to solve f_t+v f_x=0
    for i = 1:N
        for j = 1:M
            f0t[i, j] = ff0[i, j] * exp(-(2pi / L) * 1im * k_fre[j] * v[i] * t)
            f1t[i, j] = ff1[i, j] * exp(-(2pi / L) * 1im * k_fre[j] * v[i] * t)
            f2t[i, j] = ff2[i, j] * exp(-(2pi / L) * 1im * k_fre[j] * v[i] * t)
            f3t[i, j] = ff3[i, j] * exp(-(2pi / L) * 1im * k_fre[j] * v[i] * t)
        end
        f0t[i, :] .= real(ifft(f0t[i, :]))
        f1t[i, :] .= real(ifft(f1t[i, :]))
        f2t[i, :] .= real(ifft(f2t[i, :]))
        f3t[i, :] .= real(ifft(f3t[i, :]))
    end

    ################################################################
    #below we compute E1t; solve Ampere equation
    E1t = zeros(ComplexF64, M)
    E1t[1] = E1[1]
    for i = 2:(M-1)÷2+1
        #    poisson equation error in fourier version()
        #PN tmp[i] = 1im * 2pi / L * (i - 1) * E1[i] - sum(ff0[:, i]) * 2 * H / N
        #  We can also write: E1t[i] = L/(1i*2*pi*(i-1))*sum(f0t[:,i])*2*H/N
        E1t[i] =
            E1[i] +
            L / (1im * 2pi * (i - 1)) *
            sum(ff0[:, i] .* (exp.(-1im * 2pi / L * (i - 1) * t .* v) .- 1.0)) *
            2 *
            H / N
    end

    # for i = M/2+1:M
    for i = (M-1)÷2+2:M
        k = i - M - 1
        #PN tmp[i] = 1im * 2pi / L * k * E1[i] - sum(ff0[:, i]) * 2 * H / N
        E1t[i] =
            E1[i] +
            L / (1im * 2pi * k) *
            sum(ff0[:, i] .* (exp.(-1im * 2pi / L * k * t .* v) .- 1.0)) *
            2 *
            H / N
    end
    # sum(abs(tmp))
    # tmp
    return f0t, f1t, f2t, f3t, E1t
end

"""
$(SIGNATURES)
"""
function H1f!(f0, f1, f2, f3, E1, t, L, H)

    N, M = size(f0)
    #####################################################
    # computation of fO & f
    ff0 = complex(f0')
    ff1 = complex(f1')
    ff2 = complex(f2')
    ff3 = complex(f3')
    v = (1:N) .* 2 .* H ./ N .- H .- H ./ N
    # Fourier transform

    fft!(ff0, 1)
    fft!(ff1, 1)
    fft!(ff2, 1)
    fft!(ff3, 1)
    # frequency: k_fre=[0:M/2-1,-M/2:-1]
    # M is odd; odd frequency number
    k = 2pi / L .* [0:(M-1)÷2; -(M - 1)÷2:-1]
    # FFT in x direction to solve f_t+v f_x=0
    ff0 .*= exp.(-1im * k .* v' * t)
    ff1 .*= exp.(-1im * k .* v' * t)
    ff2 .*= exp.(-1im * k .* v' * t)
    ff3 .*= exp.(-1im * k .* v' * t)

    f0 .= real(ifft(ff0, 1)')
    f1 .= real(ifft(ff1, 1)')
    f2 .= real(ifft(ff2, 1)')
    f3 .= real(ifft(ff3, 1)')

    ################################################################
    #below we compute E1; solve Ampere equation
    #poisson equation error in fourier version()
    for i = 2:M
        E1[i] +=
            1.0 ./ (1im .* k[i]) .*
            sum(ff0[i, :] .* (exp.(-1im .* k[i] .* t .* v) .- 1.0)) * 2H / N
    end
end

export H1fOperator

struct H1fOperator

    adv::Translator
    tmp::Matrix{ComplexF64}

    H1fOperator(adv) = new(adv, zeros(ComplexF64, adv.mesh.M, adv.mesh.N))

end

"""
$(SIGNATURES)
"""
function step!(f0, f1, f2, f3, E1, op::H1fOperator, t)

    k = op.adv.mesh.k
    M = op.adv.mesh.M
    v = op.adv.mesh.v
    dv = op.adv.mesh.dv

    transpose!(op.tmp, f0)
    fft!(op.tmp, 1)
    op.tmp .*= exp.(-1im .* k .* v' .* t)
    for i = 2:M
        E1[i] +=
            1.0 ./ (1im .* k[i]) .*
            sum(view(op.tmp, i, :) .* (exp.(-1im .* k[i] .* t .* v) .- 1.0)) * dv
    end
    ifft!(op.tmp, 1)
    transpose!(f0, real(op.tmp))

    transpose!(op.tmp, f1)
    fft!(op.tmp, 1)
    op.tmp .*= exp.(-1im .* k .* v' .* t)
    ifft!(op.tmp, 1)
    transpose!(f1, real(op.tmp))

    transpose!(op.tmp, f2)
    fft!(op.tmp, 1)
    op.tmp .*= exp.(-1im .* k .* v' .* t)
    ifft!(op.tmp, 1)
    transpose!(f2, real(op.tmp))

    transpose!(op.tmp, f3)
    fft!(op.tmp, 1)
    op.tmp .*= exp.(-1im .* k .* v' .* t)
    ifft!(op.tmp, 1)
    transpose!(f3, real(op.tmp))
end