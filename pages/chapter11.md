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
  display_name: Julia 1.8.2
  language: julia
  name: julia-1.8
---

# Polar coordinates

```{code-cell}
include("smij-functions.jl");
```

## p28: eigenmodes of Laplacian on the disk (compare p22.jl)

```{code-cell}
# r coordinate, ranging from -1 to 1 (N must be odd):
N = 25
N2 = div(N - 1, 2)
D, r = cheb(N)
D² = D^2
D1 = D²[2:N2+1, 2:N2+1]
D² = D²[2:N2+1, N:-1:N2+2]
E1 = D[2:N2+1, 2:N2+1]
E2 = D[2:N2+1, N:-1:N2+2]

# t = theta coordinate, ranging from 0 to 2π (M must be even):
M = 20
dθ = 2π / M
θ = dθ * (1:M)
M2 = M ÷ 2
entry(k) = k==0 ? -π^2 / 3dθ^2 - 1 / 6 : 0.5 * (-1)^(k+1) / sin(dθ * k/ 2)^2
D²t = [ entry(mod(i-j,M)) for i in 1:M, j in 1:M ]

# Laplacian in polar coordinates:
R = diagm(1 ./ r[2:N2+1])
Z = zeros(M2, M2)
L = kron(D1 + R * E1, I(M)) + kron(D² + R * E2, [Z I(M2); I(M2) Z]) + kron(R^2, D²t)

# Compute four eigenmodes:
index = [1, 3, 6, 10]
λ, V = eigen(-L, sortby=abs)
λ = sqrt.(real(λ / λ[1]));
```

```{code-cell}
# Plot eigenmodes with nodal lines underneath:
using GLMakie, PyFormattedStrings

rr, θθ = r[1:N2+1], [0; θ]
XX = [r*cos(θ) for θ in θθ, r in rr]
YY = [r*sin(θ) for θ in θθ, r in rr]
z = cis.(π * (-200:200) / 200)
fig = Figure()
ax = vec([ Axis3(fig[j,i]) for i in 1:2, j in 1:2 ])
for (ax,i) in zip(ax,index)
    u = reshape(real(V[:, i]), M, N2)
    u = [zeros(M + 1) u[[M; 1:M], :]]
    u = normalize(u, Inf)
    surface!(ax, XX, YY, u)
    ax.title = f"Mode {i}:  λ = {λ[i]:.11f}"
    limits!(ax,-1.05, 1.05,-1.05, 1.05,-1.05, 1.05)
end
hidespines!.(ax)
hidedecorations!.(ax)
fig
```

### p28b: first 25 modes

can't do the contours over the disk in Makie

```{code-cell}
# Plot nodal lines:
rr, θθ = r[1:N2+1], [0; θ]
XX = [ r*cos(θ) for r in rr, θ in θθ ] 
YY = [ r*sin(θ) for r in rr, θ in θθ ] 
z = cispi.((-200:200) / 200)
fig = Figure()
ax = vec([ Axis(fig[j,i]) for i in 1:5, j in 1:5 ])
for i in 1:25
    U = reshape(real(V[:, i]), M, N2)
    U = [zeros(M + 1) U[[M; 1:M], :]]
    U = normalize(U, Inf)
    contour!(ax[i], U, levels=[0])
    ax[i].title = f"{λ[i]:.4f}"
    hidespines!(ax[i])
    hidedecorations!(ax[i])    
end
fig
```

## p29: solve Poisson equation on the unit disk

```{code-cell}
# Laplacian in polar coordinates:
N = 25   # must be odd
N2 = div(N - 1, 2)
D, r = cheb(N)
D² = D^2
D1 = D²[2:N2+1, 2:N2+1]
D2 = D²[2:N2+1, N:-1:N2+2]
E1 = D[2:N2+1, 2:N2+1]
E2 = D[2:N2+1, N:-1:N2+2]
M = 20   # must be even
dθ = 2π / M
t = dθ * (1:M)
M2 = div(M , 2)
entry(k) = k==0 ? -π^2 / 3dθ^2 - 1 / 6 : 0.5 * (-1)^(k+1) / sin(dθ * k/ 2)^2
D²t = [ entry(mod(i-j,M)) for i in 1:M, j in 1:M ]
R = diagm(1 ./ r[2:N2+1])
Z = zeros(M2, M2)
L = kron(D1 + R * E1, I(M)) + kron(D2 + R * E2, [Z I(M2); I(M2) Z]) + kron(R^2, D²t)

# Right-hand side and solution for u:
f = [ -r^2 * sin(θ / 2)^4 + sin(6θ) * cos(θ / 2)^2 for θ in θ, r in r[2:N2+1] ]
u = L \ vec(f);
```

```{code-cell}
# Reshape results onto 2D grid and plot them:
U = reshape(u, M, N2)
rr = range(0, 1, 50)
θθ = range(0, 2π, 80)
VV = zeros(80, N2)
for j in axes(U,2)
    VV[:,j] .= triginterp(U[:,j]).(θθ)
end
UU = zeros(80, 50)
for i in axes(VV,1)
    r2 = @. cos(π*(0:N)/N)
    v2 = [0; VV[i,:]; reverse(VV[mod(i-1+M2, M2) + 1,:]); 0]
    UU[i,:] .= polyinterp(r2,v2).(rr)
end

XX = [ r*cos(θ) for θ in θθ, r in rr ]
YY = [ r*sin(θ) for θ in θθ, r in rr ]
fig = Figure()
ax = Axis3(fig[1, 1], azimuth=0.5, elevation=0.55, 
    xlabel="x", ylabel="y", zlabel="z")
surf = surface!(XX, YY, UU, colormap=:bluesreds, colorrange=(-0.04,0.04))
Colorbar(fig[1,2],surf)
fig
```
