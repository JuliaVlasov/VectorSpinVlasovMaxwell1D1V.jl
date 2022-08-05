"""
computation of integral average in each cell using newton-cotes formula
$(SIGNATURES)
"""
function numeint(value, N)
    integralvalue = 7 / 90 * value[1:5:5*N-4] + 16 / 45 * value[2:5:5*N-3]
    integralvalue = integralvalue + 2 / 15 * value[3:5:5*N-2]
    integralvalue = integralvalue + 16 / 45 * value[4:5:5*N-1]
    integralvalue = integralvalue + 7 / 90 * value[5:5:5*N]

    return integralvalue
end
