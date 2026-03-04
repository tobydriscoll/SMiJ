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

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p11.jl
:label: p11
:linenos: true
:language: julia
:filename: p11
:::

### Output 11

```{code-cell}
:label: output11
using SpectralMethodsTrefethen
p11()
```

There are no new Julia tricks here, but there is a stylistic change from the original `p11.m`. In both versions, the solid curves on the left column show the given exact function. In the right column, the original program connected the dots with solid lines again, but only as a visual guide—the lines don't mean anything, really. I've had students get confused on that point, so here, I've used a dotted line to conect the points, just to show that something is different.

## Program p12

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p12.jl
:label: p12
:linenos: true
:language: julia
:filename: p12
:::

### Output 12

```{code-cell}
:label: output12
p12()
```

You're meant to compare these plots to [Output 7](#output7) for periodic functions.