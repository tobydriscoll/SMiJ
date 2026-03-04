---
title: Chapter 2
subtitle: Program p3
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
## Program p3
:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p3.jl
:label: p3
:linenos: true
:language: julia
:filename: p3
:::

### Output 3

```{code-cell}
:label: output3
using SpectralMethodsTrefethen
p3()
```

::::{note} Function definition
:icon: false
:class: dropdown
In line 7 we see one way to define a function, much like one would do it in writing, such as:

```julia
r(x) = 1 / (25x^2 + 1)
```

:::{caution}
Once a name is assigned to a function, it cannot be used as a variable name in the same scope.
:::

In lines 9–11, we see another way to define functions, with the anonymous (i.e., [lambda](wiki:Anonymous_function)) syntax

```julia
x -> float(x==0)
```
::::

::::{note} Numeric type conversion
:icon: false
:class: dropdown
In lines 9 and 10 I use `float` to convert the logical results of type `Bool` to floating-point numbers. The explicit conversions aren't strictly necessary here, but they improve clarity. A more general syntax for type conversion is `convert(T, x)`, which converts `x` to type `T`. 
::::

::::{note} `max` vs. `maximum`
:icon: false
:class: dropdown
From line 11, note that `max` chooses the largest of its arguments, unlike Chapter 1's use of `maximum`, which returns the largest element of an array.
::::

::::{note} Ternary operator
:icon: false
:class: dropdown
The {term}`ternary` operator

```julia
p ? A : B
```

evaluates to `A` if `p` is true; otherwise, `B`. This device is best used with short expressions, to keep the logic readable at a glance. Some style guides recommend against it entirely.
::::

::::{note} Subplots
:icon: false
:class: dropdown
Line 16 shows how to use indexing of a figure object to create subplots. We'll be seeing a lot of that as we go.
::::

::::{note} Generators
:icon: false
:class: dropdown
Line 18 contains more than it might seem. The goal is to evaluate a truncation of the sum in (2.10),

```{math}
\sum_{m} v_m S_h(x - x_m),
```

for all values of $x$ in the plotting range `xx`. The approach taken here is to evaluate each term at all values of $x$ and sum the resulting vectors. While the syntax resembles a {term}`comprehension`, it is instead a more primitive {term}`generator`. The distinction is clearer for a simple example:

```{code-cell}
sum( 1/factorial(n) for n in 0:8 )      # generator
```

```{code-cell}
sum( [1/factorial(n) for n in 0:8] )    # comprehension
```

These statements perform the same computation, but the former just runs a loop, while the second constructs a vector from the loop and then sums the entries of the vector.
::::

An alternative to line 18 is to sum over scalars and gather the results into a vector:

```julia
p = [sum(f(x[j]) * sinc(xx - x[j]) for j in eachindex(x)) for xx in xx]
```

It's unclear which approach is more efficient. So let's test them with the aid of `BenchmarkTools`:

```{code-cell}
using BenchmarkTools
function eval1(f) 
    v = f.(x)
    return sum( v[j] * sinc.(xx .- x[j]) for j in eachindex(x) )
end
function eval2(f) 
    return [sum(f(x[j]) * sinc(xx - x[j]) for j in eachindex(x)) for xx in xx]
end
x = -10:10
xx = -10-1/20:1/10:10+1/20
f = x -> float(abs(x) <= 3)
difference = maximum(abs, eval1(f) - eval2(f))
```

Before looking at the benchmarks, note that the two approaches are presented using the third way to define functions: with the `function` keyword. This is the preferred syntax for all but the shortest functions. The `return` keyword specifies the return value, which may be a tuple. Benchmarks are most reliable when performed by functions, rather than in the global scope. 

In order to get the best results from the `@benchmark` macro, we need to {term}`interpolate` the function `f` into the expression with `$f`. This is true for any global variable that you want to send to the benchmarking code.

```{code-cell}
@benchmark eval1($f)
```

```{code-cell}
@benchmark eval2($f)
```

From the output above, we see that the first method is significantly faster than the second, despite the number of flops being the same. When the data is small and the algorithms are fast, non-arithmetic operations can make dramatic differences.