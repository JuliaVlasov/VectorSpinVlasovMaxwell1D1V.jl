function _initialfields(x, a, ww, frequency, k0)

    E0 = 0.123 * ww # Eref
    E1 = fft(a ./ frequency .* sin.(frequency .* x))
    E2 = fft(E0 .* cos.(k0 .* x))
    E3 = fft(E0 .* sin.(k0 .* x))
    A2 = -fft(E0 ./ ww .* sin.(k0 .* x))
    A3 = fft(E0 ./ ww .* cos.(k0 .* x))

    return E1, E2, E3, A2, A3

end

function initialfields(H, L, N, M, a, ww, frequency, k0)

    x = (0:(M-1)) .* L ./ M #mesh in x direction
    _initialfields(x, a, ww, frequency, k0)

end

function initialfields(mesh, a, ww, frequency, k0)

    _initialfields(mesh.x, a, ww, frequency, k0)

end
