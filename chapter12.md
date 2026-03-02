---
title: Chapter 12
subtitle: Functions clencurt, gauss, Programs p30, p31
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p30

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p30.jl
:label: p30
:linenos: true
:language: julia
:filename: p30
:::

### Output 30

```{code-cell}
:label: output30a
using SpectralMethodsTrefethen
p30()
```

By default, [Program 30](#p30) computes quadrature weights by inverting a Chebyshev differentiation matrix, as in the text. I put that step into a function called `chebweights`, which is replaceable by the other two quadrature formulas in the chapter.

::::{note} NaN values
:icon: false
:class: dropdown
In line 26 of [Program 30](#p30), I have replaced all the zero error values with `NaN` (Not a Number). That's what causes the gaps you can see in the lower two plots of [Output 30a](#output30a). 

A `NaN` value can arise as the result of an undefined or indeterminate operation. Traditionally, it has also been used to represent missing values, and it remains convenient to use it that way in plots. But Julia also offers a separate `missing` value for that purpose, which is more appropriate for working with other kinds of data.

This brings up an amusing bit of Julia:

```{code-cell}
NaN isa Number
```

`NaN` is a special floating-point value, and all floating-point values are of type `Number`, so the assertion is true. 

:::{caution}
Numerical operations and comparisons involving `NaN` always result in `NaN`. Most notoriously,

```{code-cell}
NaN == NaN
```

To test for a `NaN` value, you should use the `isnan` function.
:::

::::


## clencurt

:::{literalinclude} SpectralMethodsTrefethen/src/clencurt.jl
:label: clencurt
:linenos: true
:language: julia
:filename: clencurt
:::

I've had some fun here by extending `clencurt` to work with floating-point types other than double precision, which is called `Float64` in Julia. Base Julia has some less-precise types, such as `Float32`, and some more-precise types, such as `BigFloat`, which is a variable-precision type. There are also add-on packages such as `DoubleFloats` and `MultiFloats` that can give higher fixed precision with greater speed than `BigFloat`. 

All that Clenshaw–Curtis requires is basic arithmetic and the cosine function, so it is easy to generalize. One has to be mindful that some operations, such as the ratio of integers, default to `Float64`, which will limit the precision of all future operations. Thus, to define `q` in line 6, I explicitly cast the numerator to have type `T`. Since the denominator is an integer, it will be implicitly converted to `T` as well before the division is performed.

Now we can integrate as accurately as we like:

```{code-cell}
T = BigFloat
setprecision(100)    # 100 bits of precision, about 30 decimal digits
println("first BigFloat after 1 is ", nextfloat(one(T)))
x, w = clencurt(40, T)
@show sum(w * exp(2x) for (x, w) in zip(x, w)) - sinh(T(2));
```

### Output 30b

```{code-cell}
:label: output30b
p30(clencurt)
```

## gauss

:::{literalinclude} SpectralMethodsTrefethen/src/gauss.jl
:label: gauss
:linenos: true
:language: julia
:filename: gauss
:::

If we want to extend `gauss` to work in higher precision, we need an eigensolver that works with that type. Those exist, but the relative cost is much higher than for Clenshaw–Curtis.

::::{note} `SymTridiagonal`
:icon: false
:class: dropdown
The `SymTridiagonal` type is a special type for symmetric tridiagonal matrices. It is used here because the matrix in the Golub–Welsch algorithm has those properties, and the eigensolver for that case is a lot more efficient than the general-purpose one.
::::

### Output 30c

```{code-cell}
:label: output30c
p30(gauss)
```

## Program p31

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p31.jl
:label: p31
:linenos: true
:language: julia
:filename: p31
:::

### Output 31

```{code-cell}
:label: output31
p31()
```
