---
title: Chapter 8
subtitle: Function chebfft, Programs p18, p19, p20
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

##  chebfft

:::{literalinclude} SpectralMethodsTrefethen/src/chebfft.jl
:label: chebfft
:filename: chebfft
:language: julia
:linenos: true
:::

I have structured my implementation of `chebfft` to take advantage of the factor of 4 efficiency mentioned in the text. The `FFTW` package gives us access to the discrete cosine and sine transforms, which I give the more convenient names `DCT` and `DST`. The DCT implicitly assumes even symmetry, which avoids the need to create an extended vector of the data. Once we differentiate cosines to get sines, we have an odd function, so the DST is the natural choice for the inverse transform; it expects us to omit the values we know to be zero at the endpoints. Everything else operates much like it does in the MATLAB code. 

::::{note} `similar`
:icon: false
:class: dropdown
The `similar` function creates an uninitialized array with the same size and element type as the given array. We have to be careful here not to specify either a real or complex array, since we want the result to have the same type as the input. We could get a valid result using `w = copy(v)`, but that causes an unnecessary copy of the data.
::::

## Program p18

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p18
using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen, ForwardDiff
"p18 - Chebyshev differentiation via FFT (compare p11)"
function p18()
    f(x) = exp(x) * sin(5x)
    fprime(x) = ForwardDiff.derivative(f, x)
    fig = Figure()
    for (i, N) in enumerate([10, 20])
        _, x = cheb(N)
        v = f.(x)

        Axis(fig[i, 1]; title=latexstring("\$u(x),\\,  N=$N\$"))
        scatter!(fig[i, 1], x, v)
        lines!(fig[i, 1], -1..1, f)

        error = chebfft(v) - fprime.(x)
        Axis(fig[i, 2]; title=latexstring("error in \$u'(x),\\,  N=$N\$"))
        scatter!(fig[i, 2], x, error)
        lines!(fig[i, 2], -1..1, chebinterp(error))
    end
    return fig
end
```

### Output 18

```{code-cell}
:label: output18
p18()
```

I made a change compared to [Output 11](#output11). For the right-hand column, I show the discrete error as markers as usual, but instead of connected them with straight lines, I use `chebinterp` to show the error as a continuous function. I like this here because it dispels the appearance that the discrete errors are local extrema.

## Program p19

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p19
using CairoMakie, LaTeXStrings, FFTW, SpectralMethodsTrefethen
"p19 - 2nd-order wave eq. on Chebyshev grid (compare p6)"
function p19(N=80, tmax=4, Δt=8/N^2)
    _, x = cheb(N)
    v = @. exp(-200x^2)
    vold = @. exp(-200 * (x - Δt)^2)

    # Time-stepping by leap frog formula:
    tplot = 0.025
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)
    data = hcat(v, zeros(N+1, ntime))
    t = Δt * (0:ntime)
    for i = 1:ntime
        w = chebfft(chebfft(v))
        w[[1, N+1]] .= 0
        vnew = 2v - vold + Δt^2 * w
        data[:, i+1] = vnew
        vold, v = v, vnew
    end

    # Plot results:
    return heatmap(x, t[1:plotgap:end], data[:, 1:plotgap:end];
        colormap=:redsblues, colorrange=(-1, 1),
        axis=(xlabel=L"x", ylabel=L"t"))
end
```

### Output 19

```{code-cell}
:label: output19
p19()
```

The reds–blues colormap makes it perfectly clear that the bump inverts whenever it reflects from a boundary. The animation is fun as well:

```{code-cell}
:tags: [remove-output]
:label: p19anim
using CairoMakie, Printf, LaTeXStrings
using FFTW, SpectralMethodsTrefethen
"p19anim - 2nd-order wave eq. on Chebyshev grid (compare p6)"
function p19anim(N=80, tmax=4)
    _, x = cheb(N)
    Δt = 8 / N^2

    tplot = 0.005
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Time-stepping by leap frog formula:
    time = Observable(0.0)
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    v = Observable(@. exp(-200x^2))
    vold = @. exp(-200 * (x - Δt)^2)
    fig = lines(x, v; axis=(; xlabel=L"x", title, limits=(-1, 1, -1, 1)))
    anim = record(fig, "p19anim.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            w = chebfft(chebfft(v[]))
            w[[1, N+1]] .= 0
            vnew = 2v[] - vold + Δt^2 * w
            vold = v[]
            v[] = vnew
            time[] = n * Δt
            iszero(mod(n, plotgap)) && recordframe!(io)
        end
      end
    return anim
end
```

```{code-cell}
:tags: [remove-output]
p19anim()
```

(output19anim)=
![](p19anim.mp4)

Compare this solution to [Output 6](#output6anim), where we had a variable speed on a periodic domain.

Your eye tells you that this solution is underresolved, but spectral methods are meant to run at coarse resolutions. To smooth out the animation, it would be easy to use the Chebyshev interpolant on a finer plotting grid.

## Program p20

```{code-cell}
:tags: [remove-output]
:class: numbered
:label: p20anim
using CairoMakie, LaTeXStrings, Printf
using FFTW, SpectralMethodsTrefethen
"p20anim - 2nd-order wave eq. in 2D via FFT (compare p19)"
function p20anim(N=24, tmax=1, Δt=6/N^2)
    # Grid and initial data:
    x = y = cheb(N)[2]
    DCT(v) = FFTW.r2r(v, FFTW.REDFT00)
    DST(v) = FFTW.r2r(v, FFTW.RODFT00)
    function fftderiv2(v)
        N = length(v) - 1
        x = [cospi(k / N) for k in 0:N]
        v̂ = DCT(v)
        dQ_dθ = DST(-(1:N-1) .* v̂[2:N]) / 2N    # omits n=0, n=N
        d²Q_dθ² = DCT(-(0:N).^2 .* v̂) / 2N
        w₂ = zero(v)    # return zero at boundaries
        # Chain rule for θ -> x:
        for n in 2:N
            s² = 1 / (1 - x[n]^2)
            s = sqrt(s²)
            w₂[n] = s² * (d²Q_dθ²[n] - x[n] * dQ_dθ[n-1] * s)
        end
        return w₂
    end

    V = [exp(-40 * ((x - 0.4)^2 + y^2)) for x in x, y in y]
    Vold = copy(V)
    Vxx = similar(V)
    Vyy = similar(V)

    tplot = 1 / 60
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Time-stepping by leap frog formula:
    time = Observable(0.0)
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    V = Observable(V)
    xx = yy = range(-1, 1, 151)
    VV = @lift interp2dgrid($V, chebinterp, chebinterp, xx, yy)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", ylabel=L"y", aspect=DataAspect(),
                limits = (-1, 1, -1, 1), title)
    heatmap!(ax, xx, yy, VV; colorrange=(-0.75, 0.75), colormap=:redsblues)
    anim = record(fig, "p20anim.mp4"; framerate=30) do io
        recordframe!(io)
        for n in 1:ntime
            Vxx = mapslices(fftderiv2, V[]; dims=[1])
            Vyy = mapslices(fftderiv2, V[]; dims=[2])
            Vnew = 2V[] - Vold + Δt^2 * (Vxx + Vyy)
            time[] = n * Δt
            Vold, V[] = V[], Vnew
            iszero(mod(n, plotgap)) && recordframe!(io)
        end
    end
    return anim
end
```

### Output 20-anim

```{code-cell}
:tags: [remove-output]
p20anim(32, 3)
```

(output20anim)=
![](p20anim.mp4)

I went straight to the animation this time. I chose to break out the FFT differentiation into its own function, `fftderiv2`. It's patterned after the `chebfft` function, but it returns the second derivative—except at the boundaries, where it returns zero. This ensures that the solution remains constant on the boundary during the time stepping. Note that the second derivative in $\theta$ is a cosine series again, so the DCT is needed to invert its transform. 

This 1-D function is then mapped to slices along each dimension in order to get the second partials needed for the wave equation. Technically, a true 2-D FFT from `FFTW` would be more efficient than repeated applications of 1-D FFTs. But the efficiency difference is unimportant unless $N$ is quite large, and the extra complication to the code is not worthwhile here.

::::{caution}
It's tempting to distill lines 28–29 into the multiple assignment

```julia
Vxx = Vyy = similar(V)
```

However, this would make both variables aliases of the same array! We need two separate arrays here.
::::

::::{note} Color range
:icon: false
:class: dropdown
The color range for the heatmap is set to $[-0.75, 0.75].$ Hence, values greater or less than these values map to the extreme colors in the colormap. While this makes the plot of the initial condition less visually informative, the amplitude of the solution diminishes quickly, and we get more color during most of the animation than we would by using the range $[-1, 1]$.
::::