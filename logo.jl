using CairoMakie, Colors
(blue, red, green, purple) = Colors.JULIA_LOGO_COLORS
f = Figure(size=(180, 90),
    fontsize = 80, fonts = (; regular = "Crimson Text Bold"),
    color = RGBAf(1, 1, 1, 0))
# ax = Axis(f[1, 1], aspect=DataAspect(), limits=((0, 3, -0.1, 1.2)))
# hidedecorations!(ax)
# hidespines!(ax)
# text!(ax, 0, 0, text="SMiJ", color = purple)
# text!(ax, 0, 0, text="SMi", color = green)
# text!(ax, 0, 0, text="SM", color = red)
# text!(ax, 0, 0, text="S", color = blue)
# resize_to_layout!(f)
chars = map(enumerate(zip("SMiJ", [blue, red, green, purple]))) do (i, (c, clr))
    rich("$c", color = clr)
end
Label(f[1, 1], rich(chars...))
f
