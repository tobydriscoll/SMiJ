---
title: Chapter 11
subtitle: Programs p28, p29
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p28

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p28.jl
:label: p28
:linenos: true
:language: julia
:filename: p28
:::

### Output 28

```{code-cell}
:label: output28a
using SpectralMethodsTrefethen
p28()
```

I found the outputs for [Program 28](#p28) in *SMiM* to be a bit coarse, so interpolated the modes to a finer grid. Here is version b:

```{code-cell}
:label: output28b
p28(version=:b)
```

Compare this result to the [modes on a square domain](#output23a).

::::{note} `mod1`
:icon: false
:class: dropdown
The most noteworthy step in setting up for interpolation on the disk is to restore the representation to all four regions in Figure 11.3 of the text. The data matrix `U` has $\theta$ along $M$ rows and $r$ along $(N-1)/2$ columns. Each row of data is for $r\in(0,1]$ only, so we have to extend it by looking in the row that is $M/2$ rows away. The trick is to treat the row index circularly, and the `mod1` function is perfect for that: it computes index modulo $M$ but with the result in $0<m\le M$ instead of $0\le m<M$. Of course, we also have to reverse the order of the columns as we append to each row. Now, the `interp2dgrid` function doesn't care that our grid is polar—it just sees a rectangular grid with our Fourier setup in the row dimension and Chebyshev in the column dimension.
::::

## Program p29