---
title: Chapter 10
subtitle: Programs p25, p26, p27
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

Here's a look at the unstable result from [Program 6](#p6):

```{code-cell}
:tags: [hide-output]
using SpectralMethodsTrefethen
p6anim(128, 6, 1.9/128)
```

![](p6anim-128-6.mp4)

You can see the same small error early on as in [the stable version](#output6anim), but that is not the real issue here. Rather, we see the high-frequency mode in space grow exponentially until it quickly takes over everything.

Here is the unstable version of p20, the [2D wave equation](#output20anim):

:::{code-cell} julia
:tags: [hide-output]
p20anim(24, 0.5, 6.6/24^2)
:::

![](p20anim-24-0.5.mp4)

You should pause the animation and step through frame by frame to see how the error sweeps in from the corners of the domain.

## Program 25

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p25.jl
:label: p25
:linenos: true
:language: julia
:filename: p25
:::

I rearranged the computations a bit for [Program 25](#p25), which shows the stability regions for various ODE solvers. For the Adams methods, I've referred to the standard characteristic polynomials $\rho$ and $\sigma$, and for the Runge–Kutta methods, I've drawn more attention to the connection to the Maclaurin series of the exponential function. 

### Output 25

```{code-cell}
:label: output25
using SpectralMethodsTrefethen
p25()
```

::::{note} `reim`
:icon: false
In the plotting commands of [Program 25](#p25), `reim(w)...` is a slight shorthand for `real(w), imag(w)`. The `reim` function returns a tuple of the real and imaginary parts of its argument, and the {term}`splatting` operator unpacks the tuple into separate arguments.

A subtle point here: the real part of an array is an array, so it's okay to use `real(w)`, although the broadcasted version `real.(w)` does the same thing. The same is true for `imag`. However, broadcasting `reim` returns a vector of tuples rather than a tuple of vectors and would not splat in a useful way here. In Makie, you could use `Point2.(reim.(w))` to get an equivalent plot. 
::::

## Program 26

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p26.jl
:label: p26
:linenos: true
:language: julia
:filename: p26
:::

The third plot is supposed to show values that oscillate around zero, but at magnitudes that are more logarithmic than they are linear:

### Output 26

```{code-cell}
:label: output26
p26()
``` 

::::{note} Flexible axis scaling
:icon: false
Makie offers [more than just linear and log scales](https://docs.makie.org/stable/reference/blocks/axis#yscale), and here I selected an inverse sinh scale:

```{math}
\tilde{y} = \frac{\mathrm{asinh} \left(y/a\right)}{\mathrm{asinh} \left(1/a\right)},
```

where $a$ is a parameter that controls the transition between linear and logarithmic behavior. The `Makie.AsinhScale` constructor takes this parameter as an argument. As you can see in [Output 26](#output26), this scale allows us to see different magnitudes of both positive and negative values.
::::

## Program 27

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p27.jl
:label: p27
:linenos: true
:language: julia
:filename: p27
:::

### Output 27

```{code-cell}
:label: output27
p27()
```

The heatmap in this case shows clearly the way the humps emerge from their interaction shifted from their original characteristic paths. That's a nonlinear effect in KdV.

I used `interpolate=true` in the `heatmap` call to get a smoother-looking plot. I guess that's a bit hypocritical given my earlier stances on keeping spectral accuracy, but I wanted to make the point that interpolating faithfully to 4th order in time is not easy here.

Simple fixed-time-step implementations of methods for IVPs are great to learn from, but they should not be used when solving challenging problems. Julia has the best collection of IVP solvers in the world, in the [DifferentialEquations](https://docs.sciml.ai/DiffEqDocs/stable/) package. Most of its sovlers automatically provide adaptive step size and interpolation in time. The package even offers [exponential Runge–Kutta methods](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#exp_RK) ready to use.

## Program 27-anim

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p27anim.jl
:label: p27anim
:linenos: true
:language: julia
:filename: p27anim
:::

For the animation, I use Fourier interpolation in space to get a smoother-looking result.

### Output 27-anim

```{code-cell}
:tags: [remove-output]
:label: output27anim
p27anim()
```

![](p27anim-256-6.0.mp4)
