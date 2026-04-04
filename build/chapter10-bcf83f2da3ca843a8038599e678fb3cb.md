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

The output shows beautifully that the soliton solutions emerge from their interaction shifted away from their otherwise straight characteristic paths in spacetime:


### Output 27

```{code-cell}
:label: output27
p27()
```

I used `interpolate=true` in the `heatmap` call to get a smoother-looking plot. That's a bit hypocritical given my earlier stances on keeping spectral accuracy, but I wanted to make the point that interpolating faithfully to 4th order in time is not easy here.

Basic fixed-time-step implementations of methods for IVPs are great for learning, but they should not be used when solving challenging problems. Julia has the best collection of IVP solvers in the world, collected in the [DifferentialEquations](https://docs.sciml.ai/DiffEqDocs/stable/) package. Most of the sovlers automatically provide adaptive step size and interpolation in time. 

## Program 27-anim

(p27anim)=
::::{code-cell} julia
:linenos: true
:filename: p27anim
using CairoMakie, LaTeXStrings, Printf
using SciMLOperators, OrdinaryDiffEq, FFTW, SpectralMethodsTrefethen
"""
p27anim - Solve KdV eq. u_t + u u_x + u_xxx = 0 on [-π, π] by
        FFT with integrating factor v = exp(-ik³t) * û.
"""
function p27anim(N=256, tmax=0.006, Δt=0.4/N^2)
    # Set up grid and two-soliton initial data:
    x = (2π / N) * (-N/2:N/2-1)
    soliton(x, a) = 3a^2 * sech(0.5 * a * x)^2
    v = @. soliton(x + 2, 25) + soliton(x + 1, 16)

    # Define linear and nonlinear parts of the PDE:
    k = [0:N/2-1; 0]
    D³ = DiagonalOperator(1im * k .^ 3)
    function nonlin(û, p, t)
        -0.5im * k .* rfft(irfft(û, N) .^ 2)
    end

    # Solve the IVP and animate solution:
    ivp = SplitODEProblem(D³, nonlin, rfft(v), (0., tmax))
    sol = solve(ivp, ETDRK4(), dt=Δt)

    time = Observable(0.0)
    title = @lift @sprintf("t = %0.2e", $time)
    v = @lift irfft(sol($time), N)
    u = @lift fourinterp($v)
    fig = lines(0..2π, u;
        axis=(; xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), title))
    anim = record(fig, "temp.mp4", range(0, tmax, 180); framerate=30) do t
        time[] = t
    end
    return anim
end
::::

For the animated form[^notpackaged] of [Program 27](#p27anim), I use the `ETDRK4` time integrator, which implements the fourth-order [exponential Runge–Kutta method](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#exp_RK). Formulating the method of lines as a `SplitODEProblem` allows us to separate the stiff linear term, which is handled implicitly, from the nonlinear term of the PDE. The resulting solution object can be called at any value of $t$ to get the solution in Fourier space.

[^notpackaged]: I did not put the `p27anim` shown here into `SpectralMethodsTrefethen`, because doing so creates some rather heavy dependencies that are otherwise not necessary. The version in `SpectralMethodsTrefethen` solves the PDE just as in [Program 27](#p27).

### Output 27-anim

```{code-cell}
:tags: [remove-output]
:label: output27anim
p27anim()
```

![](p27anim-256-6.0.mp4)

## Exercises

::::{note} Exercises 10.6, 10.7, 10.8
:icon: false
Exercises 10.6, 10.7, and 10.8 from *SMiM* can all be approached using a `SplitODEProblem`, like in [Program 27-anim](#p27anim). 
::::

::::{exercise}
:label: exr10-neumann
To repeat [Program 26](#p26) for the second derivative Chebyshev differentiation matrix with homogeneous Neumann boundary conditions, we have to do something a bit more subtle than removing the first and last rows and columns of $D_N^2$. Note that we can decompose any matrix-vector product with the discrete solution $v$ as

```{math}
\begin{aligned}
Av & = 
\begin{bmatrix}
A_{0,0} & A_{0,N} \\ 
\vdots & \vdots \\
A_{N,0} & A_{N,N}
\end{bmatrix}
\begin{bmatrix}v_0 \\ v_N
\end{bmatrix}
 + 
\begin{bmatrix}
A_{0,1} & \cdots & A_{0,N-1} \\ 
\vdots & & \vdots \\
A_{N,1} & \cdots & A_{N,N-1}
\end{bmatrix}
\begin{bmatrix}
  v_1 \\ \vdots \\ v_{N-1}
\end{bmatrix} \\ 
&= A_B v_B + A_I v_I,
\end{aligned}
```

which isolates the boundary contributions. Using the first and last rows of $D_N v= 0$ as the boundary conditions, we can rearrange and solve for $v_B = C v_I$, where $C$ is $2 \times (N-1)$. Then

```{math} 
D_N^2 v = (D_N)_B v_B + (D_N)_I v_I = [(D_N)_I + (D_N)_B C] v_I.
```

The matrix in brackets, with first and last rows dropped, is the modified differentiation matrix $\tilde{D}_N^2$ that incorporates the boundary conditions.

(a) Find the eigenvalues and eigenvectors of $\tilde{D}_N^2$ for $N=60$. Explain why the first eigenvalue (in terms of magnitude) is zero and why the second eigenvalue is close to a multiple of $\pi^2$. (You can plot their eigenvectors to get hints.)

(b) Make a plot like the top plot in [Output 26](#output26) to show the growth of the eigenvalues. Does the line at $2\pi/N$ still mark the start of a transition?

(c) Plotting an interpolated eigenvector requires first reconstructing the full vector $v$ from the interior values $v_I$. Do this to make a plot like the middle plot in [Output 26](#output26).

::::