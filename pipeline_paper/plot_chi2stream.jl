# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
"""Plot χ²Stream function in θ₀-ω₀ plane."""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using PyCall
using DelimitedFiles
using GLM
# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")
# %%


# Open input file

#chi2_file = "laruidosa/likelihood_beta0_1.258e-05.txt"
chi2_file = "chi2stream_beta0_1.254e-05.txt"
chi2_t_file = "chi2stream_tilted_beta0_1.254e-05.txt"
matriz = readdlm(chi2_file);
matriz_t = readdlm(chi2_t_file);



# %%
# Loading elements from the original run
θ = matriz[:, 1]
ω = matriz[:, 2]
β = matriz[:, 3]
χ²= matriz[:, 4] # for chi2stream file
#χ² = -matriz[:, 4]; #for likelihood file


# %%
# Looking at the Likelihood in the parameter 2D slice.
set_aog_theme!()
update_theme!(fontsize=40)
bool = χ² .< 10^4
df = (x=θ[bool], y=ω[bool], z=χ²[bool])
plt = data(df) * mapping(:x, :y; color=:z => L"χ²_{\textrm{stream}}") * visual(Scatter, markersize=1, colormap=:viridis)
bool = χ² .< 50
df_crest = (x=θ[bool], y=ω[bool], z=χ²[bool])
plt_crest = data(df_crest) * mapping(:x, :y) * visual(Scatter, markersize=5, color=:red)
sol = (36.0661, 27.3468)
fig = draw(plt + plt_crest, axis=(xlabel=L"θ_0", ylabel=L"ω_0"))
scatter!(sol, color=:black, markersize=15)
save("paper_plots/chi2stream.pdf", fig, pt_per_unit=1)
fig

# %%
set_aog_theme!()
update_theme!(fontsize=40)


f = Figure()
ax=Axis(f[1, 1], xlabel=L"θ_0", ylabel=L"ω_0", limits=((34,38),(25.5,28.2)))
levels = collect(range(1000.,10000.,length=10))
pushfirst!(levels,500)
pushfirst!(levels, 20)
co = contourf!(θ, ω, χ², levels=levels,
               extendhigh=(:black,0.2),colormap=:viridis)
Colorbar(f[1, 2], co;ticks=levels,
        ticklabelsize=25, label=L"χ²_{\textrm{stream}}")
sol = (36.033, 27.323)
scatter!(sol, color=:black, markersize=15)
save("paper_plots/chi2stream_contourf.pdf", f, pt_per_unit=1)
f

# %%
collect(range(1000.,10000.,length=10))

# %%
# Linear Fit to the valley.
linearRegressor = lm(@formula(y~x), df_crest)
a=coef(linearRegressor)
h(x)=a[1]+a[2]*x
a

# %%
# Looking at the Likelihood in the parameter 2D slice (tilted)
set_aog_theme!()
update_theme!(fontsize=40)
η = ω - h.(θ)
bool = χ² .< 1000
df_tilt = (x=θ[bool], y=η[bool], z=χ²[bool])
plt = data(df_tilt) * mapping(:x, :y; color=:z => L"χ²_{\textrm{stream}}") * visual(Scatter, markersize=10, colormap=:viridis)
bool = χ² .< 50
df_crest = (x=θ[bool], y=η[bool], z=χ²[bool])
plt_crest = data(df_crest) * mapping(:x, :y) * visual(Scatter, markersize=10, color=:red)
sol = (36.0661, 27.3468- h(36.0661))
fig = draw(plt + plt_crest, axis=(xlabel=L"θ_0", ylabel=L"ω_0-h(θ_0)",limits=((35,37),(sol[2]-0.02, sol[2]+0.02))))
scatter!(sol, color=:black, markersize=15)
save("paper_plots/chi2stream_tilted.pdf", fig, pt_per_unit=1)
display(fig)

# %%
# Loading elements from the tilted run
θ = matriz_t[:, 1]
ξ = matriz_t[:, 2]
β = matriz_t[:, 3]
χ²= matriz_t[:, 4]; # for chi2stream file


# %%
# Looking at the Likelihood in the parameter 2D slice (tilted)
# I haven't managed to make this work with AoG. See cell below with plain Makie.

set_aog_theme!()
update_theme!(fontsize=37)
df_tilt = (x=θ, y=ξ, z=χ²)
plt = data(df_tilt) * mapping(:x, :y; color=:z => L"χ²_{\textrm{stream}}") * (visual(Scatter, markersize=3, colormap=:viridis))
plt_crest = data(df_crest) * mapping(:x, :y) * visual(Scatter, markersize=1, color=:red)
sol = (36.0661, 27.3468- h(36.0661))
plt_contour =  data(df_tilt) * mapping(:x, :y, :z)*visual(Contourf,levels=[19.,19.2, 19.32, 30, 50,100,200,300,500], colormap=:viridis)

fig = Figure(fontsize = 40)
gridpos = fig[1, 1]
f = draw!(gridpos, plt_contour, axis=(xlabel=L"θ_0", ylabel=L"ω_0-h(θ_0)",limits=((35,37),(-0.01,0.01))))

AlgebraOfGraphics.colorbar!(fig[1,2],f; ticks=[19.,19.2, 19.32, 30, 50,100,200,300,500])
save("paper_plots/chi2stream_tilted.pdf", fig, pt_per_unit=1)
fig

# %%
set_aog_theme!()
update_theme!(fontsize=40)


f = Figure()
ax=Axis(f[1, 1], xlabel=L"θ_0", ylabel=L"ω_0-h(θ_0)", limits=((35,37),(-0.011,0.011)))

co = contourf!(θ, ξ, χ², levels=[19., 20, 30, 50, 70, 100,200,300,500],
               extendhigh=(:black,0.2),colormap=:viridis)

Colorbar(f[1, 2], co;ticks=[20, 50,70, 100,200,300,500],
        ticklabelsize=25, label=L"χ²_{\textrm{stream}}")
save("paper_plots/chi2stream_tilted_contourf.pdf", f, pt_per_unit=1)
f

# %%
