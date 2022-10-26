---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.1
kernelspec:
  display_name: Julia 1.8.0
  language: julia
  name: julia-1.8
---

# Differentiation matrices

## p1: convergence of fourth-order finite differences

```{code-cell} julia
using LinearAlgebra, SparseArrays

N = @. 2^(3:12)
err = zeros(size(N))
for (i,N) in enumerate(N)
    h = 2π / N
    x = @. -π + (1:N) * h
    u = @. exp(sin(x)^2)
    uʹ = @. 2 * sin(x) * cos(x) * u

    # Construct sparse fourth-order differentiation matrix:
    col1 = [ 0; -2/3h; 1/12h; zeros(N-5); -1/12h; 2/3h ]
    D = sparse( [col1[mod(i-j,N) + 1] for i in 1:N, j in 1:N] )

    # Plot max(abs(D*u-uʹ)):
    err[i] = norm(D * u - uʹ, Inf)
end
```

```{code-cell} julia
using CairoMakie

fig = Figure()
Axis(
    fig[1, 1],
    xscale=log10, yscale=log10,
    xlabel="N", ylabel="error",
    title="Convergence of fourth-order finite differences",
)

scatter!(N, err)
order4 = (N/N[1]) .^ (-4) * err[1] / 100
lines!(N, order4, color=:black, linestyle=:dash)
text!(105, 8e-8, text=L"N^{-4}")
fig
```

## p2: convergence of periodic spectral method (compare p1)

```{code-cell} julia
N = 2:2:100
@assert(all(iseven.(N)),"N must be even")
err = zeros(size(N))

for (i,N) in enumerate(N)
    h = 2π / N
    x = [-π + i * h for i = 1:N]
    u = @. exp(sin(x))
    uʹ = @. cos(x) * u

    # Construct spectral differentiation matrix:
    entry(k) = k==0 ? 0.0 : (-1)^k * 0.5cot( k * h / 2 )
    D = [ entry(mod(i-j,N)) for i in 1:N, j in 1:N ]

    # Plot max(abs(D*u - uʹ)):
    err[i] = norm(D*u - uʹ, Inf)
end
```

```{code-cell} julia
scatter(N, err, axis = (
    xlabel="N", xscale=log10, 
    ylabel="error", yscale=log10,
    title="Convergence of spectral differentiation"
))
```

