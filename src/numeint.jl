"""
computation of integral average in each cell using newton-cotes formula
$(SIGNATURES)
"""
function numeint(value, N)
    integralvalue = 7 * value[1:5:5N-4]
    integralvalue .+= 32 * value[2:5:5N-3]
    integralvalue .+= 12 * value[3:5:5N-2]
    integralvalue .+= 32 * value[4:5:5N-1]
    integralvalue .+= 7 * value[5:5:5N]

    return integralvalue ./ 90
end
