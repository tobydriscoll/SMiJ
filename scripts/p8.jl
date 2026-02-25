using LinearAlgebra, ToeplitzMatrices
"p8 - eigenvalues of harmonic oscillator -u′′ + x² u on 𝐑"
function p8(L = 8)                         # domain is [-L, L], periodic
    for N in 6:6:36
        h = 2π / N
        x = h * (1:N)
        x = L * (x .- π) / π
        column = [-π^2 / 3h^2 - 1/6; [-0.5 * (-1)^j / sin(h * j / 2)^2 for j in 1:N-1]]
        D² = (π / L)^2 * Toeplitz(column, column)    # 2nd-order differentiation
        λ = eigvals(-D² + Diagonal(x.^2))
        @show N
        println.(λ[1:4])
        println()
    end
    return nothing
end
p8()
