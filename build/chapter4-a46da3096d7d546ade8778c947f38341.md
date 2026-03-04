---
title: Chapter 4
subtitle: Programs p7, p8
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p7

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p7.jl
:label: p7
:linenos: true
:language: julia
:filename: p7
:::

### Output 7

```{code-cell}
:label: output7
using SpectralMethodsTrefethen
p7()
```

Many more new goodies appear in [Program 7](#p7).

::::{note} Range vs. vector
:icon: false
:class: dropdown

A range in Julia is not the same as a vector:

```{code-cell}
N = 6:2:50
println("N is a Vector? ", N isa Vector)
println("N is an AbstractVector? ", N isa AbstractVector)
```

Thus, while we can read an element of the range, we can't overwrite one:

```{code-cell}
:tags: [raises-exception]
@show N[4]
N[1] = 2    # error, not allowed
```

If you want to convert a range into a proper vector, you can use `collect`. It also happens implicitly in some situations, like a vertical concatenation:

```{code-cell}
println("collect(N) is a Vector? ", collect(N) isa Vector)
@show [1; N];
```

Another such situation arises in line 28. Since exponentiating the range `(-15:5:0)` creates an output that is not evenly spaced, the result is a vector.

```{code-cell}
println("10.0 ^ (-15:5:0) is a Vector? ", 10.0 .^ (-15:5:0) isa Vector)
println("10 * (-15:5:0) is a Vector? ", 10 * (-15:5:0) isa Vector)
```
::::

::::{note} `linkaxes!`
:icon: false
:class: dropdown

Line 26 uses `linkaxes!` to link the axes of all four subplots. This means that they will all have the same $x$ and $y$ limits, without the need to specify those limits explicitly.
::::

::::{note} Splatting
:icon: false
:class: dropdown
Also in line 26, we see the syntax `ax...`, which is called {term}`splatting`. This is a convenient way to indicate that all the entries of an iterable object are to be included as function arguments. ::::

::::{note} `pairs`
:icon: false
:class: dropdown

In lines 15 and 27 we have loops set up using `pairs`, like

```julia
for (i, f) in pairs(funs)
```

The `pairs` function is similar to `enumerate`, in that it returns an iterable for both the index and the value of each element of an array. However, `pairs` uses an index that is native to its argument, whereas `enumerate` always returns a range starting at 1. 

[Program 7](#p7) is set up so that `funs`, `titles`, and `ax` are all 2×2 in the same way we want them to appear in the final subplot grid. So, when we call `pairs(funs)` or `pairs(ax)`, we get an iterator whose index is a `CartesianIndex` value, like in the following:

:::{code-cell}
A = [11 12; 21 22]
for (i, v) in pairs(A)
    println("i = $i, v = $v, A[i] = $(A[i])")
end
:::

Because `i` is a `CartesianIndex`, we can use it to index into `A` directly, as in `A[i]`. Even better, we can use it as the *first two* indices into the 3-D array `E`, and then give a third index. This is incredibly convenient syntax for working with multidimensional arrays. 
::::

::::{note} Automatic differentiation
:icon: false
:class: dropdown
The `ForwardDiff` package provides a powerful way to compute derivatives called {term}`automatic differentiation`. It applies the basic rules of differentiation, like the product rule and chain rule, to compute the derivative of a Julia function at a given value of its argument. Even though it gives floating-point results, it is not a numerical approximation, but rather an exact application of the differentiation rules.

Automatic differentiation is a highly developed area for Julia, because of the benefits it offers to scientific computing (e.g., Jacobian and Hessian matrices) and machine learning (e.g., backpropagation). 
::::

::::{note} Chained assignment
:icon: false
:class: dropdown
In each of lines 30 and 31 you can see two equals signs. This is a chained assignment, which allows you to assign the same value to multiple variables in one line. They are evaluated right to left, and since the result of an assignment statement is the value that was assigned, so, for instance, the value of `ax[2, 1].ylabel` is the string `"error"`, which is then assigned to `ax[1, 1].ylabel`. 
::::

## Program p8

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p8.jl
:label: p8
:linenos: true
:language: julia
:filename: p8
:::

### Output 8

```{code-cell}
:label: output8
p8()
```

::::{note} `LinearAlgebra` and eigenvalues
:icon: false
:class: dropdown
This is our first use of the `LinearAlgebra` package, part of the standard library. While matrix multiplication and linear system solving is in base Julia, most other basic linear algebra functionality is in this package, including norms, factorizations, and eigenvalues.

Since we only want eigenvalues here, we use `eigvals` to get them. This is faster than getting both the eigenvalues and eigenvectors, which is done by `eigen`. In either case, the default is to return everything lexicographically sorted by real and imaginary parts of the eigenvalues; you can override this by adding, for example, `sortby=abs` to use magnitude (i. e., absolute value in the real case) instead.

:::{tip}
If you want to use the Unicode character λ, remember that the English spelling of this letter is "lambda", with a silent "b".
:::
::::

::::{note} Diagonal matrix
:icon: false
:class: dropdown
The `Diagonal` constructor creates a special type of matrix that stores only the diagonal elements:

```{code-cell}
using LinearAlgebra
Diagonal(1.0:5)
```

As such, it is an `AbstractMatrix`, not a `Matrix`. This allows other functions to specialize their behavior to take advantage of the diagonal structure. If you want a true `Matrix` constructed from its diagonal, use `diagm`:

```{code-cell}
diagm(1.0:5)
```

This would be necessary in order to change the off-diagonal elements, for example. Other functions will see this as a generic matrix without special structure.
::::
