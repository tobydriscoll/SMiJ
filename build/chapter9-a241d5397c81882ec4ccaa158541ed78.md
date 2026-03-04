---
title: Chapter 9
subtitle: Programs p21, p22, p23
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p21

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p21.jl
:label: p21
:linenos: true
:language: julia
:filename: p21
:::

### Output 21

```{code-cell}
:label: output21
using SpectralMethodsTrefethen
p21()
```

:::{note} `series` plot
:icon: false
:class: dropdown
The `series` function in Makie is a generalization of `lines` that can handle multiple curves at once. They correspond to the rows of the `data` matrix, which are presented in two batches, even and odd, to use different line styles.
::::

## Program p22

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p22.jl
:label: p22
:linenos: true
:language: julia
:filename: p22
:::

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

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p23.jl
:label: p23
:linenos: true
:language: julia
:filename: p23
:::

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
