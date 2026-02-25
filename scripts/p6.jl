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

    dt = h / 4
    c = @. 0.2 + sin(x - 1)^2
    v = @. exp(-100 * (x - 1) .^ 2)
    vold = @. exp(-100 * (x - 0.2dt - 1)^2)

    # Time-stepping by leap frog formula:
    tplot = 0.15
    plotgap = round(Int, tplot / dt)
    dt = tplot / plotgap
    ntime = round(Int, tmax / dt)
    data = [v zeros(N, ntime)]
    t = dt * (0:ntime)
    for i = 1:ntime
        w = fftderiv(v)
        vnew = vold - 2dt * c .* w
        data[:, i+1] = vnew
        vold, v = v, vnew
    end
    # Plot results:
    return heatmap(x, t[1:plotgap:end], data[:, 1:plotgap:end];
        colormap=:viridis,
        axis=(xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), ylabel=L"t"))
end

p6(128, 8)
