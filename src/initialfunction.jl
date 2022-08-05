
"""
Initial Gaussian function & add a perturbation in x direction

$(SIGNATURES)
"""
function initialfunction(k, x, i, v1int, frequency, a)

    kk = 0.17 # v_th
    value =
        (1 / sqrt(2 * pi) / kk) *
        exp(-(v1int[i])^2 / 2 / kk / kk) *
        (1 + a * cos(frequency * x[k]))

    return value

end
