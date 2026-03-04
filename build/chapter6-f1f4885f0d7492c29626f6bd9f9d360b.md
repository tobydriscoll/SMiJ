---
title: Chapter 6
subtitle: Function cheb, Programs p11, p12
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
(chapter6)=

## cheb

```{literalinclude} SpectralMethodsTrefethen/src/cheb.jl
:label: cheb
:filename: cheb
:language: julia
:linenos: true
```

::::{note} Short-circuit logic
:icon: false
:class: dropdown
Line 6 might look odd. It takes advantage of the {term}`short-circuiting` logical AND, `&&`. Short-circuiting means that if the first clause is false, then the second clause isn't evaluated, because the entire expression must be false anyway. This is a nice feature to avoid pesky runtime errors in conditional expressions, as in

```julia
if (length(x) >= 5) && (x[5] > 0)
```

Here, though, it's being used as a compact—some would say cryptic—alternative to the construction

```julia
if N == 0
    return zeros(1, 1), [1.0]
end
```

A similar trick works with the OR operator `||`, which skips the second clause if the first is true. For example,

```julia
(max_error < tol) || @warn("Error exceeds tolerance.")
```
::::

::::{note} `cospi`
:icon: false
:class: dropdown
The `cospi`, `sinpi`, and `tanpi` functions not only save a multiplication by $\pi$, but they can be more accurate, too, especially for large arguments.
::::

::::{note} Minimizing allocations
:icon: false
:class: dropdown
There is a subtle difference between the given line 9,

```julia
@. c[2:2:end] = -c[2:2:end]
```

and the functionally equivalent line

```julia
c[2:2:end] = -c[2:2:end]
```

The former version implies a copy of the vector from the right-hand side to the left-hand side. But the second version implies an assigment to a subarray of the vector on the left-hand side, thus requiring extra allocation of temporary storage for the subarray.

In fact, even the first version requires allocations (for the right-hand side) we don't need. The fastest version of this line is a loop:

```julia
for j in 2:2:length(c)
    c[j] = -c[j]
end
```

Coders in MATLAB and Python are (rightly) taught to avoid loops, but in Julia, loops carry no intrinsic penalty—sometimes, quite the opposite.
::::

## Program p11

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p11
using CairoMakie, LaTeXStrings
using SpectralMethodsTrefethen, ForwardDiff
"p11 - Chebyshev differentation of a smooth function"
function p11()
    f(x) = exp(x) * sin(5x)
    df_dx(x) = ForwardDiff.derivative(f, x)
    fig = Figure()
    for (i, N) in enumerate([10, 20])
        D, x = cheb(N)
        v = f.(x)

        Axis(fig[i, 1]; title=latexstring("\$u(x),\\,  N=$N\$"))
        scatter!(fig[i, 1], x, v)
        lines!(fig[i, 1], -1..1, f)

        error = D * v - df_dx.(x)
        Axis(fig[i, 2]; title=latexstring("error in \$u'(x),\\,  N=$N\$"))
        scatterlines!(fig[i, 2], x, error; linestyle=:dot)
    end
    return fig
end
```

### Output 11

```{code-cell}
:label: output11
p11()
```

There are no new Julia tricks here, but there is a stylistic change from the original `p11.m`. In both versions, the solid curves on the left column show the given exact function. In the right column, the original program connected the dots with solid lines again, but only as a visual guide—the lines don't mean anything, really. I've had students get confused on that point, so here, I've used a dotted line to conect the points, just to show that something is different.

## Program p12

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p12
using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen
"p12 - accuracy of Chebyshev spectral differentiation (compare p7)"
function p12(Nmax = 50)
    # Compute derivatives for various values of N:
    # 3rd deriv in BV, C^∞; analytic, polynomial
    funs = [x -> abs(x)^3            x -> exp(-x^(-2));
            x -> 1 / (1 + x^2)       x -> x^10]
    E = zeros(2, 2, Nmax)
    for N in 1:Nmax
        D, x = cheb(N)
        for (i, f) in pairs(funs)
            v = f.(x)
            df_dx = ForwardDiff.derivative.(f, x)
            E[i, N] = maximum(abs, D * v - df_dx)
        end
    end

    # Plot results:
    fig = Figure()
    titles = [L"|x|^3"  L"\exp(-x^2)"; L"1/(1+x^2)"  L"x^{10}"]
    ax = [Axis(fig[i, j]; title=titles[i, j], yscale=log10) for i in 1:2, j in 1:2]
    linkaxes!(ax...)
    for (i, ax) in pairs(ax)
        scatterlines!(ax, 1:Nmax, E[i, :])
    end
    ax[1, 1].ylabel = ax[2, 1].ylabel = "error"
    ax[2, 1].xlabel = ax[2, 2].xlabel = L"N"
    return fig
end
```

### Output 12

```{code-cell}
:label: output12
p12()
```

You're meant to compare these plots to [Output 7](#output7) for periodic functions.