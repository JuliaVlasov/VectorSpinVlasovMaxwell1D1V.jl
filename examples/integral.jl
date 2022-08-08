using Plots

h(v) = exp( - v^2 )

function numeint(value, N)
    integralvalue = 7 * value[1:5:5N-4]
    integralvalue .+= 32 * value[2:5:5N-3]
    integralvalue .+= 12 * value[3:5:5N-2]
    integralvalue .+= 32 * value[4:5:5N-1]
    integralvalue .+= 7 * value[5:5:5N]

    return integralvalue ./ 90
end

H = 5. / 2.
N = 129

v = (1:N) .* 2 .* H ./ N .- H
vnode = zeros(5N)

f = zeros(N)
fnode = zeros(5N)

for i = 1:N
    vnode[5i-4] = v[i] - 2H / N
    vnode[5i-3] = v[i] - (2H / N) * 3 / 4
    vnode[5i-2] = v[i] - (2H / N) / 2
    vnode[5i-1] = v[i] - (2H / N) / 4
    vnode[5i] = v[i]
end

for i = 1:5N
    fnode[i] = h(vnode[i])
end

f .= numeint(fnode, N)

plot(v, f)
plot!(vnode, fnode)
