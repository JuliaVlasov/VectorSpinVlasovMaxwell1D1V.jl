"""
$(SIGNATURES)
 
compute the first subsystem He()
- M is odd number
"""
function He(f0, f1, f2, f3, E1, E2, E3, A2, A3, t, M, N, H)

    A2t = A2 - t * E2
    A3t = A3 - t * E3
    #####################################################
    #computation of how much we should translate in the direction v1; & v2
    translatevaluev = -t .* real(ifft(E1))
    #####################################################

    for j = 1:M
        f0[:, j] .= translation((f0[:, j]), N, translatevaluev[j] .* ones(N), H)
        f1[:, j] .= translation((f1[:, j]), N, translatevaluev[j] .* ones(N), H)
        f2[:, j] .= translation((f2[:, j]), N, translatevaluev[j] .* ones(N), H)
        f3[:, j] .= translation((f3[:, j]), N, translatevaluev[j] .* ones(N), H)
    end
    
    f0t = copy(f0)
    f1t = copy(f1)
    f2t = copy(f2)
    f3t = copy(f3)

    return f0t, f1t, f2t, f3t, A2t, A3t

end

"""
$(SIGNATURES)
 
compute the first subsystem He()
- M is odd number
"""
function He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, t, H)

    A2 .= A2 .- t .* E2
    A3 .= A3 .- t .* E3
    
    e = -t .* real(ifft(E1))

    translation!(f0, e, H)
    translation!(f1, e, H)
    translation!(f2, e, H)
    translation!(f3, e, H)

end

export HeOperator

struct HeOperator

    adv :: Translator
    e :: Vector{Float64}

    function HeOperator(adv)

        e = zeros(adv.mesh.M)

        new(adv, e)

    end

end


"""
$(SIGNATURES)
 
compute the first subsystem He()
- M is odd number
"""
function step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, op :: HeOperator, t)

    A2 .-= t .* E2
    A3 .-= t .* E3
    
    op.e .= -t .* real(ifft(E1))

    translation!(f0, op.adv, op.e)
    translation!(f1, op.adv, op.e)
    translation!(f2, op.adv, op.e)
    translation!(f3, op.adv, op.e)

end
