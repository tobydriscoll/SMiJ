---
title: Chapter 5
subtitle: Programs p9, p10
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p9

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p9
using CairoMakie, Printf, Polynomials, SpectralMethodsTrefethen
"p9 - polynomial interpolation in equispaced and Chebyshev pts"
function p9(N = 16)
    xx = range(-1, 1, 801)
    fig = Figure(size=(640, 300))
    cases = [
        ("equispaced points", [-1 + 2n / N for n in 0:N]),
        ("Chebyshev points", [cospi(n / N) for n in 0:N]),
    ]
    f(x) = 1 / (1 + 16x^2)
    for (i, (title, x)) in enumerate(cases)
        v = f.(x)
        ax = Axis(fig[1, i]; title, limits=(nothing, (-1, 1.5)))
        scatter!(ax, x, v)
        u = fit(Polynomial, x, v)       # interpolation
        uu = u.(xx)                     # evaluation of interpolant
        lines!(ax, xx, uu)
        maxerr = @sprintf("%.5g", maximum(abs, uu - f.(xx)))
        text!(ax, 0, -0.5; text="max error = $maxerr", align=(:center, :bottom))
    end
    return fig
end
```

### Output 9

```{code-cell}
:label: output9
p9()
```

In *SMiM*, the notational conventions shifted a little bit for this program. Rather than using $u$ to mean a function with grid values called $v$, `p9.m` uses the variable `u` to mean grid values of a polynomial $p$. That's consistent with the math in this chapter, but not with the code in the preceding ones, and I've opted in [Program 9](#p9) to stay with $u$ to mean the function with grid values $v$. We'll see that in future programs as well.

::::{note} `range`
:icon: false
:class: dropdown
In line 4 we see a new way to create a range, with the syntax `range(begin, end, number)`. This is the analog of `linspace` in MATLAB and NumPy.
::::

::::{note} Polynomials
:icon: false
:class: dropdown

The [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) package provides ways to work with polynomials. 

:::{important}
Unlike in MATLAB, Polynomials expects coefficients to be given in order of increasing degree:

```{code-cell}
Polynomial([1, 0, -2.5, 5])
```
:::

In this case, we aren't constructing polynomials directly. Rather, I use its `fit` function to perform polynomial interpolation, analogously to `polyfit` in MATLAB. More generally, both functions do least-squares fitting if given a polynomial degree less than the lengths of the data vectors.
::::

On a subtle but important mathematical note, you would find big differences in errors between `fit` and `polyfit` at larger values of $N$. The `polyfit` function (at this writing time) computes the coefficients of the interpolating polynomial in the monomial basis, which is a bad choice due to its poor conditioning. The `fit` implementation for Polynomials, though, uses a better-conditioned approach from [](doi:10.1137/19M130100X). Of course, while that gets the interpolating polynomial accurately, the polynomial is not necessarily an accurate representation of the underlying source function when the nodes are equally spaced!

::::{note} Enumerating tuples
:icon: false
:class: dropdown
Note the syntax used with `enumerate` in line 11. Since the enumerated object is a vector whose elements are tuples, those tuples are destructured using the nested syntax `(i, (title, x))`.
::::

::::{note} Semi-automatic axis limits
:icon: false
:class: dropdown
In line 13 we see an axis being created with `limits=(nothing, (-1, 1.5))`. This syntax means to let the $x$ limits be determined automatically, while $y$ is restricted to $[-1, 1.5].$
::::

## Program p10

```{code-cell}
:class: numbered
:label: p10
using CairoMakie, Polynomials
"p10 - polynomials and corresponding equipotential curves"
function p10(N = 16)
    fig = Figure()
    cases = (
        ("equispaced points", [-1 + 2n / N for n in 0:N]),
        ("Chebyshev points", [cospi(n / N) for n in 0:N]),
    )
    for (i, (title, x)) in enumerate(cases)
        p = fromroots(x)

        # Plot p(x) over [-1,1]:
        scatter(fig[i, 1], x, zeros(size(x)), axis=(;title))
        lines!(fig[i, 1], -1..1, x -> p(x))

        # Plot equipotential curves:
        xx, yy = -1.3:0.02:1.3, -1:0.02:1
        pp = [p(x + im * y) for x in xx, y in yy]
        zz = @. log10(max(abs(pp), 1e-16))
        ax = Axis(fig[i, 2]; title, aspect=DataAspect())
        contour!(ax, xx, yy, zz; levels=-4:0, linewidth=2,
            colormap=:atlantic)
        scatter!(ax, x, zeros(size(x)); color=:black)
    end
    return fig
end
```

### Output 10

```{code-cell}
:label: output10
p10()
```

::::{note} `fromroots`
:icon: false
:class: dropdown
Line 10 uses `fromroots` to construct the monic polynomial with roots at the specified locations.
::::

::::{note} Plotting a function
:icon: false
:class: dropdown
In line 14 I introduce a tidy way to plot most functions without providing an explicit discretization of them. The syntax `lines(a..b, f)` will automatically select points in the interval $[a,b]$ to produce a smooth plot of the function $f$. This is analogous to `fplot` in MATLAB.
::::

::::{note} Avoiding log of zero
:icon: false
:class: dropdown
In line 19 I use `max` to make sure that the resulting logarithm is finite. Values of `-Inf` can confuse the contour plot.
::::

::::{note} Aspect ratio
:icon: false
:class: dropdown
In line 20, I have created an axis with specified aspect ratio. This is important in many situations, including the complex plane, where you want a unit in the $x$ and $y$ directions to have the same physical lengths in the plot. This turns out to be a [surprisingly complicated subject](https://docs.makie.org/stable/tutorials/aspect-tutorial) for plotting, depending on whether you want to adjust the axes limits (as done here with `DataAspect()`) or the shape of the axis box.
::::

::::{note} Contour plot
:icon: false
:class: dropdown
Line 20 creates the contour plot. The `levels` keyword can be an integer, to get the specified number of equally spaced levels, or a vector of explicit values. The `colormap` specifies how contour levels map to colors. You have [a lot of choices](https://docs.makie.org/stable/explanations/colors#Colormaps).
::::