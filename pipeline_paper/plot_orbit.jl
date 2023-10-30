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
"""Plot accelerations."""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using PyCall
using CSV
using DataFrames
using Showoff
using CubicSplines
using DelimitedFiles
include("bugfixMakie.jl")# %%
using Unitful, UnitfulAstro

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
galconv = pyimport("galpy.util.conversion")
galpot = pyimport("galpy.potential")
astrocoords = pyimport("astropy.coordinates")
au = pyimport("astropy.units")
importLib.reload(stream)
importLib.reload(potentials)
# importLib.reload(galpot)
# importLib.reload(galutil)
# %%

# Parameters and initial conditions.

ϵ = 56
param_file = "sol_optim_pot_m$(Int(ϵ)).txt"
θ, Δθ, β = vec(readdlm(param_file))

W = θ+Δθ
param = [ϵ, θ, W, β]


ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = vec(readdlm(ic_file))
r☼ = 8.122

orbit_nfw_file = "observable_orbit_NFW-MW.txt"
# %%

# MEPP solution
ϕ₁ₒ = stream.Iba_phi1_np
pot_list = stream.pot_model(ϵ, θ, W, β)
temp = stream.orbit_model(ic..., pot_list, r☼)
# Cubic spline for solution:
ϕ₁, ϕ₂, d☼, μ_α_cosδ, μ_δ, v☼, x, y, z, _vcirc, v_x, v_y, v_z, t = [temp[i] for i ∈ eachindex(temp)]
R = sqrt.([x[i]^2+y[i]^2 for i in eachindex(x)])
v_circ_sun = stream.rot_vel_mw(pot_list, r☼)
galcen_distance = r☼*au.kpc
v_sun = astrocoords.CartesianDifferential([11.1, v_circ_sun+12.24, 7.25]*au.km/au.s)
z_sun = 0.0*au.kpc
#%%

bool_t = t.>-0.03 .&& t.<0.03

df_b = DataFrame(ϕ₁_b=ϕ₁[bool_t],
                R_b=R[bool_t],
                z_b=z[bool_t])
df = DataFrame(ϕ₁=ϕ₁,
                R=R,
                z=z)


set_aog_theme!()
update_theme!(Axis=(topspinevisible=true, rightspinevisible=true,
topspinecolor=:darkgray, rightspinecolor=:darkgray,
xticksmirrored = true, yticksmirrored = true, aspect=1))
size_inches = (7, 7)
size_pt = 72 .* size_inches
lw = 5
fig = Figure(resolution = size_pt, fontsize = 33)
gridpos = fig[1, 1]
plt = data(df) *
    mapping(:R => L"R~[kpc]", :z => L"z~[kpc]") *
    visual(Scatter, markersize=3)
plt_b = data(df_b) *
    mapping(:R_b => L"R~[kpc]", :z_b => L"z~[kpc]") *
    visual(Scatter, markersize=6, color="red")
f = draw!(gridpos, plt+plt_b, axis=(xgridvisible=false, ygridvisible=false))
legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(55,40))
display(fig)


