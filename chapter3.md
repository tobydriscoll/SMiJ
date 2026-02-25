---
title: Chapter 3
subtitle: Programs p4, p5, p6
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## Program p4

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p4
# p4 - periodic spectral differentiation
using CairoMakie, Printf, ToeplitzMatrices
# Set up grid and differentiation matrix:
N = 24
h = 2π / N
x = h * (1:N)
col = [0.5 * (-1)^j * cot(j * h / 2) for j in 1:N-1]
D = Toeplitz([0; col], [0; reverse(col)])

fig = Figure()
xticks = MultiplesTicks(5, π, "π")

# Differentiation of a hat function:
v = @. max(0, 1 - abs(x - π) / 2)
Axis(fig[1, 1]; title="function", xticks)
scatterlines!(fig[1, 1], x, v)
Axis(fig[1, 2]; title="spectral derivative", xticks)
scatterlines!(fig[1, 2], x, D * v)

# Differentiation of exp(sin(x)):
v = @. exp(sin(x))
vprime = @. cos(x) * v
Axis(fig[2, 1]; xticks)
scatterlines!(fig[2, 1], x, v)
Axis(fig[2, 2]; xticks)
scatterlines!(fig[2, 2], x, D * v)
str = @sprintf("max error = %0.4e", maximum(abs, D * v - vprime))
text!(fig[2, 2], 2.2, 1, text=str, fontsize=12, align=(:left, :bottom))
```

### Output 4

```{code-cell}
:label: output4
fig
```

::::{note} Formatting strings
:icon: false
:class: dropdown

In [Program 4](#p4) we load a new package, `Printf`. It is part of the standard library, so it doesn't require installation. It provides the `@sprintf` macro, which enables C-style formatting of strings. You can see it in action in line 27, where `%0.4e` means to format the following number in scientific notation with four significant digits.
::::

::::{note} `reverse`
:icon: false
:class: dropdown
In line 8 you see the use of `reverse` to reverse the order of a vector. You can specify one or more dimensions for reversing an array by using the `dims` keyword argument.
::::

::::{note} `MultiplesTicks`
:icon: false
:class: dropdown
Line 11 shows off a Makie ability to create axes ticks that are multiples of any value. The call `MultiplesTicks(5, π, "π")` requests 5 ticks at multiples of $\pi$, and the tick labels will have `π` appended to them. This can be a nice refinement for plotting trigonometric functions. 
::::

::::{note} Keyword shorthand
:icon: false
:class: dropdown

Remember that arguments after the semicolon must be keyword arguments. In the calls to `Axis`, such as:

```julia
Axis(fig[1, 1]; title="function", xticks)
```

The `xticks` by itself after the semicolon is shorthand for `xticks=xticks`. When there are many keyword arguments, this trick is a nice way to shorten function calls.
::::

::::{note} `scatterlines`
:icon: false
:class: dropdown

Line 16 introduces `scatterlines!`, which combines a scatter plot (markers) with a lines plot (connecting segments).
::::

::::{note} Symbols
:icon: false
:class: dropdown

In line 28, the `align` keyword is given a tuple of two values whose names each start with a colon, `:`. This is the syntax for creating a {term}`symbol`. You can think of a symbol as a non-numeric constant, or as a string that is not meant to be printed.
::::

## Program p5

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p5
# p5 - repetition of p4 via FFT (real data only)
using CairoMakie, Printf, FFTW
# Differentiation of a hat function:
N = 24
h = 2π / N
x = h * (1:N)
function fftderiv(v)
    N = length(v)
    v̂ = rfft(v)
    ŵ = im * [0:N/2-1; 0] .* v̂
    return irfft(ŵ, N)
end

fig = Figure()
xticks = MultiplesTicks(5, π, "π")

# Differentiation of a hat function:
v = @. max(0, 1 - abs(x - π) / 2)
w = fftderiv(v)
Axis(fig[1, 1]; title="function", xticks)
scatterlines!(fig[1, 1], x, v)
Axis(fig[1, 2]; title="spectral derivative", xticks)
scatterlines!(fig[1, 2], x, w)

# Differentiation of exp(sin(x)):
v = @. exp(sin(x))
vprime = @. cos(x) * v
w = fftderiv(v)
Axis(fig[2, 1]; xticks)
scatterlines!(fig[2, 1], x, v)
Axis(fig[2, 2]; xticks)
scatterlines!(fig[2, 2], x, w)
str = @sprintf("max error = %0.4e", maximum(abs, w - vprime))
text!(fig[2, 2], 2.2, 1, text=str, fontsize=12, align=(:left, :bottom))
```

### Output 5

```{code-cell}
:label: output5
fig
```

::::{note} FFTW
:icon: false
:class: dropdown

[Program 5](#p5) is our first look at the `FFTW` package. It provides `fft` and `ifft` that are equivalent to the functions in MATLAB. But it also provides additional tools for getting greater efficiency, when doing so is important. A common case is the transformation of real data, which makes the Fourier transform an even function. The functions `rfft` and `irfft` take advantage of this symmetry to reduce the number of operations by about a factor of two. Because there is ambiguity about even/odd cases, the inverse transform needs to be given the length of the original data, as in `irfft(ŵ, N)`.
::::

::::{note} FFT differentiation
:icon: false
:class: dropdown

Unlike in the original book, I've chosen to put the operation of differentiation via FFT into a short function. Until recently, it was impossible to define a multiline function at an arbitrary point within a MATLAB script, so it was easier to simply repeat short motifs as needed. But it's better practice to use functions for this: it clarifies the purpose of the code, makes it more reusable, and eliminates bugs due to copy–paste errors or customizations. In Julia, this structure can also improve performance, for technical reasons.
::::

::::{note} Hat accent
:icon: false
:class: dropdown

In math, including *SMiM*, it's conventional to put a hat over the name of a function to indicate its Fourier transform. In Julia, you can type {kbd}`v\hat` {kbd}`TAB` to get `v̂`, for example. Note that the hat is just a character, and it doesn't have any special meaning to Julia.

Such accents can be hard to read in some fonts, so this practice comes down to personal preference.
::::

::::{note} Imaginary unit
:icon: false
:class: dropdown

In line 11 we see a reference to `im`, the imaginary unit. In Julia, `im` is a built-in constant with the value $\sqrt{-1}$. You can also write `1im` or `2im`, for example, to get $i$ or $2i$. Notice that just as Julia draws distinctions between integers and floating-point numbers, it also distinguishes between real and complex numbers:

```{code-cell}
println("Type of 1 is ", typeof(1))
println("Type of 1.0 is ", typeof(1.0))
println("Type of 1im is ", typeof(1im))
println("Type of 1.0im is ", typeof(1.0im))
```
::::

::::{note} Mixing numeric types
:icon: false
:class: dropdown

Operations involving mixed numeric types try to promote the operands to a common type. This often goes without notice. But attempting a type incorrectly can produce errors:

```{code-cell}
:tags: [raises-exception]
x = [1, 2, 3]
x[2.0]    # error, can't use floating-point number as an index
```

So can trying to store a value into an array of a less-general type:

```{code-cell}
:tags: [raises-exception]
z = [π, 2π, 3π]
@. z *= 1im    # error, can't store complex in a real array
```

:::{tip}
When you see an `InexactError`, it may be the result of a numeric type that can't be used in context without changing its value.
:::
::::

## Program p6

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p6
using CairoMakie, LaTeXStrings, FFTW
"p6 - variable coefficient wave equation"
function p6(N, tmax)
    # Grid, variable coefficient, and initial data:
    h = 2π / N
    x = h * (1:N)
    function fftderiv(v)
        N = length(v)
        v̂ = rfft(v)
        ŵ = im * [0:N/2-1; 0] .* v̂
        return irfft(ŵ, N)
    end

    Δt = h / 4
    c = @. 0.2 + sin(x - 1)^2
    v = @. exp(-100 * (x - 1) .^ 2)
    vold = @. exp(-100 * (x - 0.2Δt - 1)^2)

    # Time-stepping by leap frog formula:
    tplot = 0.15
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)
    data = [v zeros(N, ntime)]
    t = Δt * (0:ntime)
    for i = 1:ntime
        w = fftderiv(v)
        vnew = vold - 2Δt * c .* w
        data[:, i+1] = vnew
        vold, v = v, vnew
    end
    # Plot results:
    return heatmap(x, t[1:plotgap:end], data[:, 1:plotgap:end];
        colormap=:viridis,
        axis=(xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), ylabel=L"t"))
end
```

### Output 6

```{code-cell}
:label: output6
p6(128, 8)
```

Our first PDE solution!

::::{note} `heatmap`
:icon: false
:class: dropdown

You will certainly have noticed that the output of the program is not at all like Output 6 in *SMiM*. The `waterfall` plot in MATLAB is not easy to do compactly in Makie. Rather than try to reconstruct it, I have opted to use a `heatmap` instead, which conveys the result equally well, if a bit less lusciously.
::::

::::{note} Function wrapper
:icon: false
:class: dropdown

You probably wonder why I wrapped all the action into a function. The proximate answer is complicated as well as boring, and I go into it at [the end of this chapter](#why-function). But the broader answer is that Julia is oriented towards working with functions most of the time. Functions don't interfere with each others' values—except through mutation, which is supposed to be made transparent through the bang convention—or with global state. Julia's compiler is best at optimizing functions, too.

The preference for functions is so significant that:

:::{important}
 From this point onward, all the programs are written as functions.
::: 

The biggest downside to working this way is with debugging, which is one of Julia's weak points. At this writing, my recommendation is to use the [Infiltrator](https://github.com/JuliaDebug/Infiltrator.jl) package, which allows you to copy variables from a function's workspace into the global Main workspace for later inspection.
::::

::::{note} Documentation string
:icon: false
:class: dropdown

Line 2 demonstrates a documentation string for a function. A string value immediately preceding a function keyword provides a string that is printed out at the prompt when you enter `?` followed by the function name.
::::

::::{note} `round`
:icon: false
:class: dropdown

Note in lines 21 and 23 that the default output type of `round` is the same as the input type, so to use the result as an array index, we need to convert it to an integer with `round(Int, ...)`.
::::

::::{note} Horitontal concatenation
:icon: false
:class: dropdown

In line 24, I use square brackets to construct an array, with horizontal concatenation implied by the space. In Julia, vectors are one-dimensional arrays, but when a vector is used in a matrix context, it is treated as a column vector. Thus, the result of `[v zeros(N, ntime)]` is an array whose first column is `v` and whose remaining `ntime` columns are zeros. You could also use `hcat` or `cat(.., dims=2)` for horizontal concatenation.
::::

## Animated alternative to Output 6

One advantage I have here over a printed book is the use of animations. Accordingly, here is an alternative version of [Program 6](#p6) that produces an animation of the solution as it evolves in time.

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p6anim
using CairoMakie, LaTeXStrings, FFTW
"p6anim - variable coefficient wave equation"
function p6anim(N, tmax)
    h = 2π / N
    x = h * (1:N)
    function fftderiv(v)
        N = length(v)
        v̂ = rfft(v)
        ŵ = im * [0:N/2-1; 0] .* v̂
        return irfft(ŵ, N)
    end

    Δt = h / 4
    c = @. 0.2 + sin(x - 1)^2

    # Time-stepping by leap frog formula:
    ntime = ceil(Int, tmax / Δt)
    Δt = tmax / ntime
    time = Observable(0.0)
    v = Observable(@. exp(-100 * (x - 1) .^ 2))
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    vold = @. exp(-100 * (x - 0.2Δt - 1)^2)
    fig = lines(x, v;
        axis=(; xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), title))
    anim = record(fig, "p6anim.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            w = fftderiv(v[])
            vnew = vold - 2Δt * c .* w
            vold = v[]
            time[] = n * Δt
            v[] = vnew
            recordframe!(io)
        end
    end
    return anim
end
```

### Output 6-anim

```{code-cell}
:tags: [remove-output]
p6anim(128, 12);
```

(output6anim)=
![](p6anim.mp4)

Before digging into the code, look at what the original [Output 6](#output6) was hiding! While the solution in *SMiM* looks virtually perfect, we see in the animation (as well as the heatmap, if you look closely) that the solution has small but significant flaws. The culprit is the leapfrog time-stepping: at second order, it's much less accurate than the Fourier differentiation in space. Moreover, its numerical error propagates without dissipation. Many other time-stepping methods introduce some artificial dissipation, which would make the solution look a lot better—but not necessarily more accurate.

::::{note} Animation with Makie
:icon: false
:class: dropdown
The animation mechanism in Makie is to define so-called Observable quantities whose changing values drive the animation. In [Program 6-anim](#p6anim), there are two Observable values: scalar `time` and the discrete solution `v`. Other variables whose values are derived from an Observable quantity should be defined using the `@lift` macro, as in line 21 to define `title`. Plot objects created using Observable values, or values lifted from them, will be updated whenever the underlying Observable is changed. Thus, our `lines` object in lines 23–24 depends both on `v` and on `title`, which in turn depends on `time`.

The animation is created by a call to `record`. The only remaining quirk is that when you want to get or change the value of an Observable, you have to "dereference" it by adding empty square brackets to the name. You can see I have done this to every use of `v` and `time` within the recording loop.
::::

::::{note} `latexstring`
:icon: false
:class: dropdown
One more fine point in this code is the use of `latexstring` in line 21. This is another function from the `LaTeXStrings` package. We can't use the `L"` string macro here, because we need to LaTeX-ify a string whose value is unknown at compile time due to the interpolation of the time value. Note that the `\$` is needed in the string to get a literal dollar sign for LaTeX to process.
::::

::::{note} Output suppression
:icon: false
:class: dropdown
In the call to `p6anim`, I have added a semicolon to suppress the output. If you use MATLAB, you know that any line that doesn't end with a semicolon will print the value of the expression. That's true when running Julia interactively—but only for the last value of a code block. Here, `p6anim` just returns a string with the name of the file that was created. I've suppressed that and then had the web page display the resulting animation.
::::

(why-function)=
## Why use a function?

Remember that loops have locally scoped variables, except for variables that are already in the outer scope at the beginning of the loop.  The use of local variable scopes in loops (and other blocks) is rather uncommon. In other languages, including MATLAB, a reference to a variable inside a loop always refers to a value that is the same before, within, and after the loop. But in Julia, the behavior inside a loop depends on what has come before. In [Program 6](#p6), lines 16–17 define values for `v` and `vold`, so inside the time-stepping loop, those names should refer to the outer variables, *not* new local ones.

However, for historical reasons, these rules change in exactly one case: at the global scope of an interactive session. If you were to copy the interior of the function and paste it into an interactive session, you would have to add the statement

```julia
global v, vold
```

within the loop to get the same behavior. This declaration woud not cause any problems to have in any case, but the code is just cleaner if we put everything into a function scope and avoid [the whole dumb mess](https://docs.julialang.org/en/v1/manual/variables-and-scoping/).
