---
title: Chapter 7
subtitle: Programs p13, p14, p15, p16, p17
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p13

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p13
using CairoMakie, Printf, LinearAlgebra, SpectralMethodsTrefethen
"p13 - solve linear BVP u_xx = exp(4x), u(-1) = u(1) = 0"
function p13(N = 16)
    D, x = cheb(N)
    D² = D^2
    D² = D²[2:N, 2:N]                    # boundary conditions
    rhs = @. exp(4x[2:N])
    v = D² \ rhs                         # Poisson eq. solved here
    v = [0; v; 0]
    xx = range(-1, 1, 801)
    u = chebinterp(v)                    # interpolate grid data
    exact(x) = @. (exp(4x) - sinh(4) * x - cosh(4)) / 16;
    maxerr = maximum(abs, u.(xx) - exact.(xx))
    title = @sprintf("max err = %.3g", maxerr)
    fig = lines(-1..1, u; axis=(; title, xlabel=L"x", ylabel=L"u(x)"))
    scatter!(x, v)
    fig
end
```

### Output 13

```{code-cell}
:label: output13
p13()
```

::::{note} Superscript
:icon: false
:class: dropdown
Why not use a superscript for the second derivative matrix? You get it by typing {kbd}`D\^2`, then {kbd}`TAB`.
::::

I have departed from the original `p13.m` in a significant way. Rather than using `polyfit` to interpolate the solution to a finer grid, I have used the barycentric interpolation formula presented in the exercises of Chapter 6. This is a much faster and stabler method, as Trefethen made clear in a [paper published after *SMiM*](doi:10.1137/S0036144502417715).

Here is the implementation of general barycentric interpolation:

:::{literalinclude} SpectralMethodsTrefethen/src/baryinterp.jl
:label: baryinterp
:filename: baryinterp
:language: julia
:linenos: true
:::

The implementation is split into computation of the weights and application of the formula. This allows the caller to provide precomputed weights if they are available, which is the case for Chebyshev points. The rescaling in line 8 of `baryweights` is discussed in the paper I cited above.

::::{note} Conditional generator
:icon: false
:class: dropdown
Line 9 shows that a generator can include an `if` statement, which makes the entire weight generation process a one-liner.
::::

::::{note} Optional arguments
:icon: false
:class: dropdown
The definition of `baryinterp` in line 20 shows how to define an optional input argument. If only two arguments are given, then `w` is computed by a call to `baryweights`. All earlier arguments are available at this step.
::::

::::{note} Returning a function
:icon: false
:class: dropdown
The return value of `baryinterp` is a function, the interpolant. The syntax in line 23 sets this up as an anonymous function. It will have access to all the variables in scope at the time of its defintion, including the enclosing scope.
::::

::::{note} `nothing`
:icon: false
:class: dropdown
The `nothing` in line 24 is a special value of type `Nothing`. Nothing equals `nothing` except `nothing`, though the syntax `isnothing` is more idiomatic to check for it. I'm using it here to keep track of whether the given value of `s` is equal to one of the nodes. Mathematically, that creates a L'Hôpital's rule situation, but the code has to check for it explicitly.
::::

::::{note} `break`
:icon: false
:class: dropdown
The `break` statement in line 29 is used to exit the loop immediately when a node is found that matches the input `s`. There's no need to finish the loop, since the value of the interpolant is known at this point.
::::

::::{note} `dot`
:icon: false
:class: dropdown
The `dot` function in line 32 is defined in the `LinearAlgebra` package, and it computes the dot product of two vectors. In the complex-valued case, it does the proper inner product by taking the complex conjugate of the first vector. 

This is a good opportunity to point out a consequence of the fact that all vectors in Julia have a column orientation. It might surprise you that the following snippet fails:

```{code-cell}
:tags: [raises-exception]
A = ones(3, 3)
A[1, :] * A[:, 3]    # error, both operands are 3×1
```

Instead, you could use a conjugate transpose, `A[1, :]' * A[:, 3]`, or you could just use `dot`.
::::

With `baryinterp` in hand, we can specialize it to the Chebyshev case, where we don't even need to provide the nodes:

:::{literalinclude} SpectralMethodsTrefethen/src/chebinterp.jl
:label: chebinterp
:filename: chebinterp
:language: julia
:linenos: true
:::

## Program p14

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p14
using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen
"""
p14 - solve nonlinear BVP u_xx = exp(u), u(-1) = u(1) = 0 (compare p13)
"""
function p14(N = 16)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    v = zeros(N-1)
    it = 0
    while true              # fixed-point iteration
        vnew = D² \ exp.(v)
        change = maximum(abs, vnew - v)
        v = vnew
        it += 1
        (change < 1e-15 || it > 99) && break
    end
    v = [0; v; 0]
    p = chebinterp(v)
    fig = Figure()
    title = "no. steps = $it,  u(0) = $(p(0))"
    ax = Axis(fig[1, 1]; title, xlabel=L"x", ylabel=L"u(x)")
    scatter!(ax, x, v)
    lines!(ax, -1..1, p)
    return fig
end
```

### Output 14

```{code-cell}
:label: output14
p14()
```

I use `chebinterp` to get the value of the solution at $x=0$, rather than indexing into the discrete solution directly. This would work just as well for odd values of $N$, where there is no node at $x=0$.

## Program p15

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p15
using CairoMakie, Printf
using SpectralMethodsTrefethen, LinearAlgebra
"p15 - solve eigenvalue BVP u_xx = λ u, u(-1) = u(1) = 0"
function p15(N = 36)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    λ, V = eigen(-D²)
    fig = Figure(size=(600, 600))
    ax = [Axis(fig[i, 1], limits=(-1, 1, -0.64, 0.64)) for i in 1:6]
    for (mode, ax) in zip(5:5:30, ax)        # plot 6 eigenvectors
        v = [0; V[:, mode]; 0]
        scatter!(ax, x, v)
        u = chebinterp(v)
        lines!(ax, -1..1, u)
        eig = "eig $mode = $(-λ[mode] * 4/π^2) π²/4"
        text!(ax, 0, 0.45; text=eig, fontsize=12, align=(:center, :baseline))
        PPW = @sprintf("%.1f", 4N / (π * mode))
        text!(ax, 0.7, 0.45; text="$PPW  ppw", fontsize=12, align=(:left, :baseline))
    end
    hidespines!.(ax)
    hidedecorations!.(ax)
    return fig
end
```

### Output 15

```{code-cell}
:label: output15
p15()
```

::::{note} Eigenvalues and eigenvectors
:icon: false
:class: dropdown
Since we want both eigenvalues and eigenvectors, we use the `eigen` function rather than `eigvals`.
::::

::::{note} Figure size
:icon: false
:class: dropdown
The `size` keyword in line 9 of `p15` allows us to specify the size of the figure. The exact meaning of the units is [a bit complicated](https://docs.makie.org/stable/explanations/figure#Figure-size-and-resolution), but my goal here is to get more vertical space than the default of 600 by 450.
::::

::::{note} `zip`
:icon: false
:class: dropdown
The `zip` function in line 10 allows us to loop over two identically sized collections at once.
::::

::::{note} `hidespines!` and `hidedecorations!`
:icon: false
:class: dropdown
The `hidespines!` and `hidedecorations!` functions in line 17 are used to hide the axes and ticks, which makes the plot look cleaner and matches the original style. Note that they are vectorized to be applied to all the axes at once, though we could have called them from within the loop instead.
::::

## Program p16

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p16
using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"p16 - Poisson eq. on [-1,1] x [-1,1] with u = 0 on boundary"
function p16(N = 24)
    # Set up grids and tensor product Laplacian and solve for u:
    D, x = cheb(N)
    y = x
    f(x, y) = 10sin(8x * (y - 1))
    F = [f(x, y) for x in x[2:N], y in y[2:N]]
    D² = (D^2)[2:N, 2:N]
    Δ = kron(I(N-1), D²) + kron(D², I(N-1))    # Laplacian

    fig = Figure(size=(640, 800))
    ax1 = Axis(fig[1, 1]; title="Nonzeros in the Laplacian", 
      aspect=1, yreversed=true)
    spy!(ax1, Δ)
    @time v = Δ \ vec(F)    # solve problem and watch the clock

    # Reshape long 1D results onto 2D grid (reversing to usual direction):
    V = zeros(N+1, N+1)
    V[2:N, 2:N] .= reshape(v, N-1, N-1)
    icen = div(3N, 4) + 1
    value = V[icen, icen]

    # Interpolate to finer grid and plot:
    xx = yy = range(-1, 1, 201)
    UU = interp2dgrid(V, chebinterp, chebinterp, xx, yy)
    M = maximum(abs, UU)
    ax2 = Axis3(fig[2, 1]; title="Solution of Poisson equation",
      xlabel=L"x", ylabel=L"y", zlabel=L"u(x,y)")
    surface!(ax2, xx, yy, UU; colormap=:redsblues, colorrange=(-M, M))
    str = latexstring("u(2^{-1/2},2^{-1/2}) = $(@sprintf("%0.11f", value))")
    text!(ax2, -1, 1, 1.4M; text=str, align=(:left, :top))
    rowsize!(fig.layout, 1, Relative(0.35))
    return fig
end
```

### Output 16

```{code-cell}
:label: output16
p16()
```

A lot is going on here! This is our first look at a problem with two space variables. My comments on the interpolation aspect of [Program 16](#p16) are given [in detail below](#interp2d), but first, here are notes on the rest of the program.

::::{note} Nested loop in a comprehension
:icon: false
:class: dropdown
The syntax in line 9 is comprehension with a nested loop in a comprehension. There's no need to repeat the `for` keyword. Like most things in a comprehension, the ordering is backwards from writing things out longhand: 

```julia
for (j, y) in enumerate(y[2:N])
    for (i, x) in enumerate(x[2:N])
        F[i, j] = f(x, y)
    end
end
```

Thus, the *first* loop in the comprehension varies the most quickly.

:::{warning}
History unfolded such that in a matrix, the first dimension is vertically oriented, while for a bivariate function, the first dimension is horizontally oriented. Because MATLAB predates widespread use of graphics-cabable displays, it adopted the convention that in a matrix that represents the discretization of a bivariate function, the first dimension is vertically oriented. Matplotlib, being patterned after MATLAB, adopted the same convention. In modern contexts, this choice can create a lot of confusion, and it's particularly unfortuitous when there are more than two independent variables.

In Julia, the `Plots` package follows the MATLAB convention, while `Makie` more naturally has the first dimension (rows) map to the first variable ($x$). If you see codes that swap the placement of $x$ and $y$ from what I use, or that have transposed matrices in plotting calls, that's probably what's going on.
:::
::::

::::{note} Identity matrix
:icon: false
:class: dropdown
The `I` in line 11 is defined in the `LinearAlgebra` package. Given on its own, `I` represents the identity matrix of whatever size is needed for the context. Here, where the size cannot be inferred, I have given it explicitly. Keep in mind that `I` is an abstract matrix, not a true one; in some situations, you may want to convert it, as in `Matrix(I(n))`.
::::

::::{note} Kronecker product
:icon: false
:class: dropdown
The `kron` function, defined in `LinearAlgebra`, computes the Kronecker product of two matrices. My use of row dimension for $x$ and column dimension for $y$ means that `kron(I(N-1), D²)` corresponds to the second derivative in $x$, while `kron(D², I(N-1))` corresponds to the second derivative in $y$. Since we treat the dimensions identically here, of course, the distinction doesn't actually matter.
::::

::::{note} Δ
:icon: false
:class: dropdown
If being able to use Δ as the name of the discrete Laplacian doesn't convince you to use Unicode in Julia, then I don't know what to tell you. You can get it by typing {kbd}`\Delta`, then {kbd}`TAB`. Better yet, define a text replacement rule or a keyboard mapping for it.
::::

::::{note} Spy plot
:icon: false
:class: dropdown
The `spy!` function in line 13 creates a plot in Makie to show thesparsity pattern of a matrix. Note that I supplied `yreversed=true` in the axis creation to make the vertical dimension increase going downward, like row index does. 
::::

::::{note} Timing
:icon: false
:class: dropdown
The `@time` macro in line 17 reports execution time and memory allocation for the following expression. It's part of base Julia, which is convenient. But it's not as reliably reproducible as `BenchmarkTools` if you need comparative timing information.
::::

::::{note} Reshaping
:icon: false
:class: dropdown
We have two inverse types of reshaping in this program, to go between a matrix representing values on a grid and its interpretation as a linear vector. We could do both using `reshape`, but I have used `vec` to go from matrix to vector.

Be advised that reshaping does not copy data on its own, so you can end up with two variables that point to the same data, as in this example:

```{code-cell}
a = ones(4 * 6 * 10)
A = reshape(a, :, 6, 10)    # one dimension size can be left out
a[31] = 0
findfirst(iszero, A)
```

There is a more efficient but clunkier alternative to line 21:

```julia
vv = reshape(v, N-1, N-1)
copyto!(V, CartesianIndices((2:N, 2:N)), vv, CartesianIndices((1:N-1, 1:N-1)));
```
::::

::::{note} Integer division
:icon: false
:class: dropdown
The `div` function in line 22 performs integer division (quotient without remainder). You can write it in infix form as `3N ÷ 4`, typing {kbd}`\div` then {kbd}`TAB` to get the symbol. To get the remainder, you can use `rem`, and to get both at once, you can use `divrem`.
::::

::::{note} Surface plot
:icon: false
:class: dropdown
A surface is drawn in an `Axis3` object. For plots of functions that go roughly the same above and below $z=0$, I like to use a divergent colormap, like `:redsblues`, and to set the color range to be symmetric about zero, as in `colorrange=(-M, M)`.
::::

::::{note} `rowsize!`
:icon: false
:class: dropdown
If you don't indicate otherwise, each subplot in a grid layout gets the same amount of space. In line 34 I use `rowsize!` to give the first row just 35% of the total height of the figure.
::::

(interp2d)=
### Interpolation in 2-D

I have chosen to use 2-D Chebyshev interpolation to get data on the plotting grid. You might not notice the difference from a low-order method, such as bicubic spline interpolation, on a plot, but the numerical values would be different—and the low-order ones would be much less accurate. Fortunately, it is not hard to turn one-dimensional interpolation into 2-D interpolation on a tensor product grid, as I show in the `interp2dgrid` function from `SpectralMethodsTrefethen`. 

Suppose we have an interpolation method $F$ that turns values on a grid $x_1,\ldots,x_m$ into a function $f$,

```{math}
f(\xi) = F(v_1, \ldots, v_m; \xi),
```
and another method $G$ that does the same for a grid $y_1,\ldots,y_n$,

```{math}
g(\eta) = G(w_1, \ldots, w_n; \eta).
```

Then with a matrix of grid values, $v_{ij} = u(x_i, y_j)$, we can define a 2-D interpolant by a nested application of the 1-D methods:

```{math}
\begin{aligned}
u(\xi, \eta) &= G(w_1(\xi), \ldots, w_n(\xi); \eta),\\ 
 w_j(\xi) &= F(v_{1j}, \ldots, v_{mj}; \xi), \quad j = 1, \ldots, n.
\end{aligned}
```
Evaluation on a grid of values $\xi_1, \ldots, \xi_p$ and $\eta_1, \ldots, \eta_q$ leads to

```{math}
\begin{aligned}
W_{ij} = w_j(\xi_i) &= F(v_{1j}, \ldots, v_{mj}; \xi_i), \quad i = 1, \ldots, p,\\
U_{ij} = u(\xi_i, \eta_j) &= G(W_{i1}, \ldots, W_{in}; \eta_j), \quad j = 1, \ldots, q.
\end{aligned}
```

Here is my compact implementation of this idea:

```{literalinclude} SpectralMethodsTrefethen/src/interp2dgrid.jl
:label: interp2dgrid
:filename: interp2dgrid
:language: julia
:linenos: true
```

The first call to `mapslices` applies the first interpolation method to each column of data, evaluating at the given values of $\xi$, to get the matrix $W$ defined above. The second call applies the second interpolation method to each row of $W$, evaluating at the given values of $\eta$, to get the final matrix $U$.

## Program p17

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p17
using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"""
p17 - Helmholtz eq. u_xx + u_yy + k²u = f on [-1,1] x [-1,1] (compare p16)
"""
function p17(N = 24, k = 9)
    # Set up spectral grid and tensor product Helmholtz operator:
    D,x = cheb(N)
    y = x
    F = [ exp(-10*((y - 1)^2 + (x - 0.5)^2)) for x in x[2:N], y in y[2:N] ]
    D² = (D^2)[2:N,2:N]
    L = kron(I(N-1), D²) + kron(D², I(N-1)) + k^2 * I    # Helmholtz

    # Solve for u, reshape to 2D grid, and plot:
    v = L \ vec(F)
    V = zeros(N+1, N+1)
    V[2:N, 2:N] = reshape(v, N-1, N-1)
    xx = yy = range(-1, 1, 201)
    UU = interp2dgrid(V, chebinterp, chebinterp, xx, yy)

    fig = Figure(size=(600, 800))
    ax1 = Axis3(fig[1, 1]; title="Solution of Helmholtz equation",
      xlabel=L"x", ylabel=L"y", zlabel=L"u(x,y)")
    M = maximum(abs, UU)
    surface!(ax1, xx, yy, UU; colormap=:redsblues, colorrange=(-M, M))
    icen = div(N, 2) + 1
    text!(ax1, -1, 1, M; text="u(0, 0) = $(@sprintf("%0.11f", V[icen, icen]))")
    ax2 = Axis(fig[2, 1]; xlabel=L"x", ylabel=L"y", aspect=1)
    contour!(ax2, xx, yy, UU; levels=range(-M, M, 11), colormap=:redsblues)
    lines!(ax2, [-1, -1, 1, 1, -1], [-1, 1, 1, -1, -1], color=:black, linewidth=2)
    hidedecorations!(ax2; label=false, ticks=false)
    return fig
end
```

### Output 17

```{code-cell}
:label: output17
p17()
```

Nothing new needed here; just enjoy the view.
