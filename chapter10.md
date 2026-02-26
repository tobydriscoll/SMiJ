---
title: Chapter 10
subtitle: Programs p25, p26, p27
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

Here's a look at the unstable version of p6:

```{code-cell}
:tags: [hide-output]
using SpectralMethodsTrefethen
p6anim(128, 6, 1.9/128)
```

![](p6anim-128-6.mp4)

You can see the same small error early on as in [the stable version](#output6anim), but that is not the real issue here. Rather, we see the high-frequency mode in space grow exponentially until it quickly takes over everything.

Here is the unstable version of p20, the [2D wave equation](#output20anim):

:::{code-cell}
:tags: [hide-output]
p20anim(24, 0.5, 6.6/24^2)
:::

![](p20anim-24-0.5.mp4)

You should pause the animation and step through frame by frame to see how the error creeps in from the corners of the domain.

## Program p25

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p25.jl
:label: p25
:linenos: true
:language: julia
:filename: p25
:::

### Output 25

```{code-cell}
:label: output25
using SpectralMethodsTrefethen
p25()
```

I rearranged the computations a bit. For the Adams methods, I've referred to the standard characteristic polynomials $\rho$ and $\sigma$, and for the Runge–Kutta methods, I've drawn more attention to the connection to the Maclaurin series of the exponential function. 

::::{note} `reim`
I've used `reim(w)...` as a slight shorthand for `real(w), imag(w)`. The `reim` function returns a tuple of the real and imaginary parts of its argument, and the {term}`splatting` operator unpacks the tuple into separate arguments.

A subtle point: the real part of a vector is a vector, so it's okay to use `real(w)`, though the broadcasted version `real.(w)` would also work. The same is true for `imag`. However, broadcasting `reim` returns a vector of tuples rather than a tuple of vectors and would not splat in a useful way here. In Makie, you could use `Point2.(reim.(w))` to get an equivalent plot. 
::::

## Program p26

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p26.jl
:label: p26
:linenos: true
:language: julia
:filename: p26
:::

### Output 26

```{code-cell}
:label: output26
p26()
``` 

I didn't connect the markers in the third plot above because the eigenmode's behavior is a great deal more complicated than that suggests: it changes sign frequently, which would mean zooming to negative infinity and back on the log scale in the plot. 

## Program p27

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

## Program p27-anim

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p27anim.jl
:label: p27anim
:linenos: true
:language: julia
:filename: p27anim
:::

### Output 27-anim

```{code-cell}
:label: output27anim
p27anim()
```

![](p27anim-256-6.0.mp4)
