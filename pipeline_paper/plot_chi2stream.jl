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
#     display_name: Julia 1.9.0
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
# chi2_file = "chi2stream_beta0_1.258e-05.txt"
chi2_file = "laruidosa/likelihood_beta0_1.258e-05.txt"
matriz = readdlm(chi2_file);




# %%
θ = matriz[:, 1]
ω = matriz[:, 2]
β = matriz[:, 3]
#χ²= matriz[:,4] # for chi2stream file
χ² = -matriz[:, 4]; #for likelihood file


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
# Looking at the Likelihood in the parameter 2D slice (tilted)

set_aog_theme!()
size_inches = (6.2*2, 4*2)
size_pt = 72 .* size_inches
fig = Figure(resolution = size_pt, fontsize = 40)

η = ω - h.(θ)
bool = χ² .< 1000
df_tilt = (x=θ[bool], y=η[bool], z=χ²[bool])
plt = data(df_tilt) * mapping(:x, :y; color=:z => L"χ²_{\textrm{stream}}") * (visual(Scatter, markersize=10, colormap=:viridis))
bool = 40 .< χ² .< 50
df_crest = (x=θ[bool], y=η[bool], z=χ²[bool])
plt_crest = data(df_crest) * mapping(:x, :y) * visual(Scatter, markersize=10, color=:red)
sol = (36.0661, 27.3468- h(36.0661))
fig = Figure(resolution = size_pt, fontsize = 37)
gridpos = fig[1, 1]
f = draw!(gridpos,plt+plt_crest, axis=(xlabel=L"θ_0", ylabel=L"ω_0-h(θ_0)",limits=((35,37),(sol[2]-0.02, sol[2]+0.02))))
scatter!(sol, color=:black, markersize=25)

save("paper_plots/chi2stream_tilted.pdf", fig, pt_per_unit=1)
display(fig)

# %%
