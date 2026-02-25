using CairoMakie, Printf
using SpectralMethodsTrefethen, LinearAlgebra
"p15 - solve eigenvalue BVP u_xx = λ u, u(-1) = u(1) = 0"
function p15(N = 36)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    λ, V = eigen(-D²)
    fig = Figure()
    ax = [Axis(fig[i, 1], limits=(-1, 1, -0.64, 0.64)) for i in 1:6]
    for (j, ax) in zip(5:5:30, ax)        # plot 6 eigenvectors
        v = [0; V[:, j]; 0]
        scatter!(ax, x, v; markersize=5)
        u = chebinterp(v)
        lines!(ax, -1..1, u)
        eig = "eig $j = $(-λ[j] * 4/π^2) π²/4"
        text!(ax, 0, 0.5; text=eig, fontsize=8, align=(:center, :baseline))
        PPW = @sprintf("%.1f", 4N / (π * j))
        text!(ax, 0.7, 0.5; text="$PPW  ppw", fontsize=8, align=(:left, :baseline))
    end
    hidespines!.(ax)
    hidedecorations!.(ax)
    return fig
end
p15()
