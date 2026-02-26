---
title: Chapter 8
subtitle: Function chebfft, Programs p18, p19, p20
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

## chebfft

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

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p18.jl
:label: p18
:linenos: true
:language: julia
:filename: p18
:::

### Output 18

```{code-cell}
:label: output18
using SpectralMethodsTrefethen
p18()
```

I made a change compared to [Output 11](#output11). For the right-hand column, I show the discrete error as markers as usual, but instead of connected them with straight lines, I use `chebinterp` to show the error as a continuous function. I like this here because it dispels the appearance that the discrete errors are local extrema.

## Program p19

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p19.jl
:label: p19
:linenos: true
:language: julia
:filename: p19
:::

### Output 19

```{code-cell}
:label: output19
p19()
```

The reds–blues colormap makes it perfectly clear that the bump inverts whenever it reflects from a boundary. The animation is fun as well:

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p19anim.jl
:label: p19anim
:linenos: true
:language: julia
:filename: p19anim
:::

```{code-cell}
:tags: [remove-output]
p19anim()
```

(output19anim)=
![](p19anim-80-4.mp4)

Compare this solution to [Output 6](#output6anim), where we had a variable speed on a periodic domain.

Your eye tells you that this solution is underresolved, but spectral methods are meant to run at coarse resolutions. To smooth out the animation, it would be easy to use the Chebyshev interpolant on a finer plotting grid.

## Program p20-anim

:::{literalinclude} SpectralMethodsTrefethen/src/scripts/p20anim.jl
:label: p20anim
:linenos: true
:language: julia
:filename: p20anim
:::

### Output 20-anim

```{code-cell}
:tags: [remove-output]
p20anim(32, 3)
```

(output20anim)=
![](p20anim-32-3.mp4)

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