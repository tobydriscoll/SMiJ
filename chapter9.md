---
title: Chapter 9
subtitle: Programs p21, p22, p23
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p21

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p21
using CairoMakie, LaTeXStrings, LinearAlgebra, ToeplitzMatrices
"""
p21 - eigenvalues of Mathieu operator -u_xx + 2 q cos(2x) u\n
         (compare p8 and p. 724 of Abramowitz & Stegun)
"""
function p21(N = 42)
    h = 2π / N
    x = h * (1:N)
    c0 = -π^2 / 3h^2 - 1 / 6
    col = [c0; [ (-1)^(j+1) / 2sin(h * j / 2)^2 for j in 1:N-1 ]]
    D² = Toeplitz(col, col)    # second-order differentiation matrix
    q = 0:.2:15
    data = zeros(11, length(q))
    for (i, q) in enumerate(q)
        λ = eigvals( -D² + 2q * Diagonal(cos.(2x)) )
        data[:, i] = λ[1:11]
    end
    fig = Figure(size=(320, 500))
    ax = Axis(fig[1, 1]; xlabel=L"q", ylabel=L"λ", yticks=-24:4:32,
        limits=(0, 15, -24, 32))
    series!(ax, q, data[1:2:end, :], color=:atlantic)
    series!(ax, q, data[2:2:end, :], linestyle=:dot, color=:atlantic)
    return fig
end
```

### Output 21

```{code-cell}
:label: output21
p21()
```

:::{note} `series` plot
:icon: false
:class: dropdown
The `series` function in Makie is a generalization of `lines` that can handle multiple curves at once. They correspond to the rows of the `data` matrix, which are presented in two batches, even and odd, to use different line styles.
::::

## Program p22

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p22
using CairoMakie, LaTeXStrings, Printf
using LinearAlgebra, ToeplitzMatrices, SpecialFunctions, SpectralMethodsTrefethen
"p22 - 5th eigenvector of Airy equation u_xx = λ*x*u"
function p22()
    fig = Figure(size=(600, 600))
    N = [12 24; 36 48]
    for (i, N) in pairs(N)
        D, x = cheb(N)
        D² = (D^2)[2:N, 2:N]
        λ, V = eigen(D², diagm(x[2:N]))         # generalized ev problem
        fifth = findfirst(>(0), λ) + 4
        λ, v = λ[fifth], V[:, fifth]
        v = [0; v; 0]                           # extend to booundary
        v = v / v[div(N, 2) + 1] * airyai(0)    # normalize
        u = chebinterp(v)
        lines(fig[i[1], i[2]], -1..1, u;
            axis=(; title=@sprintf("N = %d,  λ = %0.10f", N, λ)))
    end
    return fig
end
```

### Output 22

```{code-cell}
:label: output22
p22()
```

::::{note} Generalized eigenvalue problem
:icon: false
:class: dropdown
The `eigen` function can solve the generalized eigenvalue problem $Av = \lambda Bv$. If $B$ is nonsingular, you can rewrite this as an ordinary eigenvalue problem $B^{-1}Av = \lambda v$, but it is not always desirable to do so. However, that is the algorithm used if $B$ is a `Diagonal` matrix, so here we have to use a generic matrix instead to avoid an error due to the singularity of $B$.
::::

::::{note} Special functions
:icon: false
:class: dropdown
The [SpecialFunctions](https://specialfunctions.juliamath.org/stable/functions_overview/) package provides a large collection of special functions, such as the gamma function, Bessel functions, elliptic integrals, the error function, the zeta function, and, as used here, the Airy functions.
::::

:::{note} `findfirst`
:icon: false
:class: dropdown
The `findfirst` function returns the index of the first element of an array that satisfies a given condition. The result is either `nothing`, if no elements satisfy the condition, or an integer index. I have used the syntax `>(0)` to create an anonymous function that tests whether its argument is greater than zero, shorthand for `x -> x > 0`. Finally, we want the 5th index of the result, so I add 4 since they are in sorted order.
::::

## Program p23

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p23
using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"""
# p23 - eigenvalues of perturbed Laplacian on [-1,1] x [-1,1] (compare p16)
"""
function p23(N = 16; perturb=false)
    # Set up tensor product Laplacian and compute 4 eigenmodes:
    D, x = cheb(N)
    y = x

    D² = (D^2)[2:N, 2:N]
    Δ = -kron(I(N-1), D²) - kron(D², I(N-1))    # Laplacian
    if perturb
        F = [exp(20 * (y - x - 1)) for x in x[2:N], y in y[2:N]]
        Δ += Diagonal(vec(F))
    end
    λ, V = eigen(Δ)

    # Reshape them to 2D grid, interpolate to finer grid, and plot:
    xx = yy = range(-1, 1, 101)
    U = zeros(N+1, N+1)
    fig = Figure(size=(600, 600))
    modes = [1 2; 3 4]
    ax = [Axis(fig[i, j]; aspect=1) for i in 1:2, j in 1:2]
    for (i, mode) in pairs(modes)
        U[2:N, 2:N] = reshape(V[:, mode], N-1, N-1)
        UU = interp2dgrid(U, chebinterp, chebinterp, xx, yy)
        UU /= maximum(abs, UU)
        str = @sprintf("λ = %0.12f π^2/4", 4λ[mode] / π^2)
        contourf!(ax[i], xx, yy, UU; levels=-1.1:0.2:1.1, colormap=:redsblues)
        ax[i].title = str
        hidedecorations!(ax[i], ticks=false)
    end
    return fig
end
```

### Output 23a

Compare the outputs to [Output 16](#output16), which was a solution of a boundary-value problem rather than an eigenvalue problem. 

```{code-cell}
:label: output23a
p23()
```

### Output 23b

```{code-cell}
:label: output23b
p23(perturb=true)
```

Both versions of `p23` are so similar that I have put them into a single function with a keyword argument to control the perturbation.

::::{note} `contourf`
:icon: false
:class: dropdown
The `contourf` function in Makie is used to create filled contour plots. The result is a kind of hybrid between a heatmap and a contour plot.
::::

::::{note} Axis title
:icon: false
:class: dropdown
Rather than adding the title at the time of axis creation, I have added it by assigning to the `.title` property of the axis after the fact.
::::
