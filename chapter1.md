---
title: Chapter 1
subtitle: Programs p1, p2
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p1

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p1
# p1 - convergence of fourth-order finite differences
using CairoMakie, LaTeXStrings
using ToeplitzMatrices
# For various N, set up grid in [-π, π] and function u(x):
N = [2^j for j in 3:12]
errors = []
for N in N
    h = 2π / N
    x = [-π + j * h for j in 1:N]
    u = @. exp(sin(x)^2)
    uprime = @. 2sin(x) * cos(x) * u

    # Construct sparse fourth-order differentiation matrix:
    col = zeros(N)
    col[[2, 3, N-1, N]] = [-2/3, 1/12, -1/12, 2/3] / h
    D = Toeplitz(col, -col)

    # Compute max of abs(D*u - uprime):
    error = maximum(abs, D * u - uprime)
    push!(errors, error)
end

fig = Figure()
ax = Axis(fig[1, 1]; title="Convergence of fourth-order finite differences",
    xscale=log10, yscale=log10, xlabel=L"N", ylabel="error")
scatter!(ax, N, errors)
lines!(ax, N, 1 ./ N .^ 4; linestyle=:dash, color=:gray)
text!(ax, 105, 1e-8; text=L"N^{-4}")
```

### Output 1

```{code-cell}
:label: output1
fig
```

There is so much to cover as we get started! These explanations are wordy, because I'm trying to be transparent, but you don't necessarily have to absorb everything at once in order to understand what's happening. I've summarized the main takeaway at the beginning of each note.

::::{note} Importing packages
:class: dropdown
:icon: false

:::{tip}
The `using` keyword imports a package. So does `import`, with different implications on how things are referred to by name.
:::

Let's begin with the opening lines,

```julia
using CairoMakie, LaTeXStrings
using ToeplitzMatrices
```

Only a small core called Base is loaded into Julia at startup. To get more functionality, we need to load packages with `using` or `import`.  The first line, for example, loads the `ToeplitzMatrices` package, provided it has been installed. Each package creates its own {term}`namespace` that holds all the objects it defines. Most packages also choose to {term}`export` some objects into the default `Main` workspace. Here, this includes the `Toeplitz` function that is called later. But even if it had not been exported, it could be accessed as the fully qualified name `ToeplitzMatrices.Toeplitz`. A package needs to be loaded only once per session. For clarity, I will be explicitly loading all packages needed at least once per chapter, so that you can see the dependencies.

The `CairoMakie` package refers to one of the two most popular packages for making plots, known as [Makie](https://docs.makie.org/stable/). (The other popular choice is [Plots](https://docs.juliaplots.org/latest/).) Makie has a few different backends for writing output in different contexts; CairoMakie is the most "print-friendly" variant. Finally, `LaTeXStrings` makes it easy to construct a string that will be interpreted as LaTeX when used in a plot.
::::

::::{note} Comprehensions, ranges, and vectors
:class: dropdown
:icon: false

```{important}
Square brackets creates create arrays from what's inside them. Ranges have the same colon syntax as MATLAB.
```

Now we move to line 5,

```julia
N = [2^j for j in 3:12]
```

This syntax is known as a {term}`comprehension`, which will be familiar to Python users. Essentially, it's a compact way of writing a loop, with the body of the loop given first. The enclosing square brackets indicate that the result will be an array. The expression `2^j` is evaluated for each value of `j` in the range `3:12`, resulting in a one-dimensional array with entries $2^3, 2^4, \ldots, 2^{12}$.

:::{#def-range}
The syntax `a:b` creates a range of values from `a` to `b`, inclusive of both endpoints. The values have a step size of 1 by default, but this can be changed with a third argument, as in `a:step:b`.
:::

Arrays are a major part of Julia. The array type is called `Array`, and it is parameterized by the type of its entries and the number of dimensions:

```{code-cell}
typeof(N)
```

As you see here, a `Vector` is simply a one-dimensional array.

Line 6 creates an empty array called `errors` that will be used to accumulate approximation errors in the loop.
::::

::::{note} Loops
:class: dropdown
:icon: false

```{caution}
Loops have a local variable scope that can overlap with the enclosing scope. This can take some getting used to.
```

The opening of the loop in line 7 looks curious:

```julia
for N in N
```

The loop variable `N` is local to the loop and distinct from the `N` in its enclosing scope. It will be assigned each value within the outer `N` in turn, and then vanish when the loop ends. This idiom is a way to avoid having to have a different name such as `N_values` for the array. If you don't like this convention, you are free to use any two different names here, in which case the outer variable will be accessible inside the loop.

This aspect of Julia can cause confusion. For example, the following code snippet, on its own, causes an error:

```{code-cell}
:tags: [raises-exception]
for i in 1:5
    n = 2i
end
println("n = ", n)    # error, n is not defined
```

Since `n` went out of scope at the end of the loop, it is undefined when the print statement executes. You can get around this by declaring `n` to be `global` within the loop, or by giving `n` a value before the loop starts, in which case the loop variable now refers to the outer `n`.
::::

::::{note} Syntactic sugar
:class: dropdown
:icon: false

:::{tip}
Greek letters, superscripts and subscripts, many symbols, and even emoji can be parts of variable names. You can also write `4x` to mean `4 * x`.
:::

At line 7 we see

```julia
h = 2π / N
```

The symbol `π` is a built-in constant in Julia. This is a good place to point out one cosmetic but really nice feature of Julia: it allows Unicode characters. You can type `π` with the shortcut `\pi<TAB>`, and it will be treated as a single character. (Personally, I use Greek letters and other symbols so often that I have many alt-key combinations mapped to them.) If you can't be bothered to type `\pi<TAB>`, you can also write `pi`, but that isn't true in general. For instance, `σ` is not the same as `sigma`.

Also notice above that we write `2π` without a multiplication sign, just as you would mathematically. You see this implied multiplication in line 10 as well, with `2sin(x)`. This won't work for multiplying two variables together, though, because, e.g., `xy` is a valid variable name, so you would need to use `x * y`.

Julia was written with mathematics in mind from the start, and it is full of small delights like these.
::::

::::{note} Vectorization
:class: dropdown
:icon: false

```{important}
Use a period to make any function be applied elementwise to an array. Use `@.` to do this to everything within the expression.
```

In line 9, we see a {term}`macro`, which is signaled by a name starting with the `@` character. A macro is a piece of code that runs at compile time and transforms the code that follows it. Here, the `@.` macro is a convenient way to {term}`vectorize` an expression.

Mathematical functions such as `sin` and addition are intended to be applied to scalar values. In scientific computing, though, we often find ourselves wanting to apply a function to every element of a vector or other array. In MATLAB and NumPy, most predefined mathematical functions will accept an array as input and do exactly that. Julia does *not*, so the expressions `sin([0, π/2, π])` and `1 + [1, 2, 3]` give errors. This can be frustrating, though such errors may expose bugs.

But there is a huge upside. You can put a period `.` at the end of a function name (or in front of an infix operator like `+`) to have it be interpreted elementwise. Unlike MATLAB and NumPy, this works for *any* function, including one you define.

It gets better. Suppose you vectorize an expression to get something like

```julia
sin.(x).^2 .- 1
```

The Julia parser knows to "fuse" these vectorizations together into a single loop, rather than creating a vector for `sin.(x)`, another vector for squaring, and a third for subtracting 1. This is great for performance, but it's not so easy on the eyes. That is what `@.` is for; it implies vectorization wherever needed:

```julia
@. sin(x)^2 - 1
```

The simplicity and speed of vectorization is one of Julia's greatest strengths for scientific computing.[^exp]

[^exp]: There is a mathematical angle here, too: sometimes, functions are defined differently for scalars and other objects. The classic example is the exponential function. In Julia, `exp(A)` computes the matrix exponential (not elementwise) if `A` is a square matrix, and throws an error if it's rectangular.
::::

::::{note} Array indexing
:class: dropdown
:icon: false

```{important}
Indexing starts at 1.[^begin] The last index is the keyword `end`.
```

Line 13 constructs a vector of $N$ zeros. Then, in line 14, we get

```julia
col[[2, 3, N-1, N]] = [-2/3, 1/12, -1/12, 2/3] / h
```

This is a way to assign values to specific entries of the vector `col`. The left-hand side refers to elements at indices 2, 3, $N-1$, and $N$, and the right-hand side is a 4-vector as well.

[^begin]: While starting at index 1 is the behavior of standard arrays, there are packages that create array-like objects that can be indexed more arbitrarily. The absolute safest way to get the first index in general is to use the `begin` keyword.
::::

::::{note} Abstract types
:class: dropdown
:icon: false

:::{tip}
Things can be array-like without being literal arrays in memory. Most of the time, you won't be aware of the difference. 
:::

Line 16 creates the finite-difference matrix `D` by calling the `Toeplitz` function. Every row in a Toeplitz matrix is a circularly shifted version of the one above it, and so the matrix is fully specified by its first column and first row:

```{code-cell}
col = [0; -2; 6; zeros(3); -6; 2]
A = Toeplitz(col, -col)
```

Perhaps surprisingly, the object `A` is *not* a matrix:

```{code-cell}
A isa Matrix
```

However, it behaves like an array in most ways. For instance, we can get its size:

```{code-cell}
size(A)
```

In line 18 of `p1`, we also see that matrix-vector multiplication is defined for a `Toeplitz` matrix with the `*` operator. Julia has a special way of indicating that an object is matrix-like, even when it is not represented by a contiguous block of memory:

```{code-cell}
A isa AbstractMatrix
```

We won't need to refer to this concept often, but it's worth understanding, because Julia relies deeply on its type system.
::::

::::{note} Minimizing allocations
:class: dropdown
:icon: false

:::{tip}
When you have a need for speed, Julia gives you access to fine-grained control.
:::

Looking more closely at line 19, we see that the error is computed by taking the difference of two vectors, `D * u` and `uprime`, and then applying the `maximum` function to its entries while applying `abs`, the absolute value, as a filter. We would get the same result from

```julia
maximum(abs.(D * u - uprime)) 
```

However, the syntax in [Program 1](#p1) is slightly preferable, because it avoids the need to create an intermediate vector of absolute values.

Julia often lets you achieve things in multiple ways. Often, it doesn't matter which you choose. But if a key section of code is taking more time than you like for manipulations like memory management that are not inherently needed, Julia gives you ways to try alternatives.
::::

::::{note} Mutating functions
:class: dropdown
:icon: false

```{caution} 
A function with a bang `!` at the end of its name can be expected to alter the values or states of its input arguments.
```

In line 20, we use `push!` to append a new value to the end of the `errors` array. The `!` at the end of a function name is a convention that indicates that the function is {term}`mutating` its arguments.

Technically, this is a little wasteful, since we have to allocate a new vector each time through the loop. We'll look at a tidier version in the next program.
::::

::::{note} Plotting
:class: dropdown
:icon: false

```{important}
Keyword function arguments come after positional ones.
```

The rest of the program creates the plot. We use `Figure` to create a figure andthen  `Axis` to create a 2D axis set within it:

```julia
ax = Axis(fig[1, 1]; title="Convergence of fourth-order finite differences",
    xscale=log10, yscale=log10, xlabel=L"N", ylabel="error")
```

We'll see later that indexing the figure object lets us create grid layouts. After the first argument, the semicolon indicates that the remaining arguments are given by keyword rather than position. (You can often use a comma instead, but the semicolon is never wrong.) Here, the keywords describe aspects of how the axis object is to look.

The actual data plotting is done by `scatter!` (unconnected markers) and `lines!` (unmarked positions connected by line segments), which both mutate the axis object. Line 28 creates a text object that will be formatted by LaTeX thanks to the `L` before the opening quote.[^Lmacro]

[^Lmacro]: This is a shortcut for invoking a macro called `@L_str`. You may see other string macros too, like `r"` for regular expressions.

Having `fig` as the result of a line causes the figure to be displayed. This has to be done every time you make changes to a figure.
::::

## Program p2

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p2
# p2 - convergence of periodic spectral method (compare p1)
using CairoMakie, LaTeXStrings, ToeplitzMatrices
# For various N (even), set up grid as before:
N = 2:2:100
errors = zeros(size(N))
for (k, N) in enumerate(N)
    h = 2π / N
    x = [-π + j * h for j in 1:N]
    u = @. exp(sin(x))
    uprime = @. cos(x) * u

    # Construct spectral differentiation matrix:
    col = [0.5 * (-1)^i * cot(i * h / 2) for i in 1:N-1]
    D = Toeplitz([0; col], [0; reverse(col)])

    # Compute max(abs(D*u - uprime)):
    error = maximum(abs, D * u - uprime)
    errors[k] = error
end

fig = Figure()
ax = Axis(fig[1, 1]; title="Convergence of spectral differentiation",
    xscale=log10, yscale=log10, xlabel=L"N", ylabel="error")
scatter!(ax, N, errors)
```

### Output 2

```{code-cell}
:label: output2
fig
```

::::{note} enumerate
:class: dropdown
:icon: false

:::{tip}
Use `enumerate` to get both index and value as you loop though a vector. Use parentheses to create a tuple, which is immutable.
:::

In [Program 2](#p2), which is structurally similar to [Program 1](#p1), we avoid unnecessary intermediate allocations for the error vector by allocating all the necessary space before the loop begins. Now, however, we need to keep track of the index `k` where the next error value should be stored. The `enumerate` function is a convenient way to get both the index and the value of each entry in the loop. The syntax `(k, N)` indicates a {term}`tuple` of two values. By putting the tuple on the left-hand side of an assignment, we can {term}`destructure` the tuple on the right-hand side, as in the following:

```{code-cell}
a, b, c = (1, 2, 3)
a + b + c
```

The main functional difference between a vector and a tuple is that a vector is mutable, so its entries can be changed, while a tuple is immutable.
::::

::::{note} Array concatenation
:class: dropdown
:icon: false

```{important}
Row dimension is the first dimension. Separate columns with spaces and rows with semicolons.
```

The syntax `[0; col]` in line 11 is worth a mention. Inside square brackets, commas separate entries of a vector. Because the elements of a vector can be anything, including other vectors, this implies that we cannot use commas to splice scalars and vectors into a single vector the way we do with scalars alone:

```{code-cell}
col = [1, 2, 3]    # a vector of integers
[0, col]           # a vector of two entries: 0 and col
```

Instead, we use a semicolon to get a vertical concatenation of the scalar and the vector. Note that the following are equivalent:

```julia
[0; col]
vcat(0, col)
cat(0, col, dims=1)    # row dimension is the first dimension
```
::::