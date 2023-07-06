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
"""Rotation curves computation and plotting."""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using PyCall
using DelimitedFiles
using CSV
using DataFrames
# %%

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
galpy_pots = pyimport("galpy.potential")
u = pyimport("astropy.units")
importLib.reload(stream)

# %%
# Parameters and initial conditions.
param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
θ, ω, β = readdlm(param_file)

ϵ = 56.0
W = θ + ω
param = [ϵ, θ, W, β]


ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122

orbit_nfw_file = "observable_orbit_NFW-MW.txt"
# %%
# Fermionic-MW solution and rotation curve
pot_list = stream.pot_model(ϵ, θ, W, β)
halo = pot_list[4]
r = 10 .^ range(log10(halo.r_s[begin]),log10(100.),length=2000)
v_f = [stream.rot_vel_mw(pot_list, x) for x in r]

# NFW-MW solution and rotation curve
function malhan_vcirc(r)
    bp= galpy_pots.PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.0,normalize=0.05)
    mp= galpy_pots.MiyamotoNagaiPotential(a=3.0/8.0,b=0.28/8.0,normalize=0.6)
    hp= galpy_pots.TriaxialNFWPotential(a=16.0/8.0,b=1.0,c=0.82,normalize=0.59)
    mw = bp+mp+hp
    v0=mw[1].vcirc(r*u.kpc)
    v1=mw[2].vcirc(r*u.kpc)
    v2=mw[3].vcirc(r*u.kpc)
    return sqrt(v0*v0+v1*v1+v2*v2)*220.0
end
v_nfw = [malhan_vcirc(x) for x in r];
df_models = DataFrame(r=r, v_nfw=v_nfw, v_f=v_f)

# %%
# Reading rotation curve observables
# Set rotation data

df_Sof13 = DataFrame(CSV.File("data_rotcurves/vel_Sofue13.txt"; delim=" ", ignorerepeated=true))
df_Sof20 = DataFrame(CSV.File("data_rotcurves/vel_Sofue20.txt"; delim=" ", ignorerepeated=true))
df_Eilers = DataFrame(CSV.File("data_rotcurves/vel_Eilers.txt"; delim=" ", ignorerepeated=true))
# asym_error = DataFrame([v_Eilers.e_down, v_Eilers.e_up])
df_Sof13.r = df_Sof13.r / 1.0e3
df_Sof13.err_r = df_Sof13.err_r / 1.0e3
df_Sof20

# %%

# Plot
labels = ["Fermionic-MW", "NFW-MW"]
lw = 4

#let
set_aog_theme!()
size_inches = (6.2 * 2, 3 * 2)
size_pt = 72 .* size_inches
fig = Figure(resolution=size_pt, fontsize=37)
gridpos = fig[1, 1]
grp = dims(1) => renamer(labels) => ""
grp2 = dims(1) => renamer(["Sofue"]) => ""
plt_model = data(df_models) *
        mapping(:r => L"r~[\textrm{kpc}]", [2, 3] .=> L"v_{\textrm{circ}}~[\textrm{km s}^{-1}]";
            color=grp,
            linestyle=grp
        ) *
        visual(Lines, linewidth=lw)
plt_Sof13 = data(df_Sof13) * mapping(:r, :v) * visual(Scatter, direction=:x)
plt_eb_Sof13 = data(df_Sof13) * (mapping(:r, :v, :err_r) * visual(Errorbars, direction=:x) + mapping(:r, :v, :err_v) * visual(Errorbars))
plt_eb_Sof20 = data(df_Sof20) * mapping(:Radius, :Velocity, :Error) * visual(Errorbars)
plt_full = plt_model+plt_Sof13+plt_eb_Sof13+plt_eb_Sof20
f = draw!(gridpos, plt_full, axis=(; limits=((0, 40), (0, 300)),
    xgridvisible=false, ygridvisible=false))
legend!(gridpos, f; tellwidth=false, halign=:center, valign=:bottom, margin=(10, 10, 10, 10), patchsize=(50, 35))
# Lines re-styling
amber_aog = "#ffa700"
green_aog = "#107A78"
lineas = fig.content[1].scene.plots
lineas[1].linestyle = :dash
lineas[2].linestyle = :dot
lineas[1].color = amber_aog
lineas[2].color = green_aog
lineas[3].color = :blue

leg = fig.content[2]
_lines = leg.blockscene.children[1].plots[2:3]
for l in _lines
    l.linewidth = 4
end
_lines[1].linestyle = :dash
_lines[2].linestyle = :dot
_lines[1].color = amber_aog
_lines[2].color = green_aog

display(fig)
save("paper_plots/rotation_curves.pdf", fig, pt_per_unit=1)
println("plot done.")
#end

# %%
leg.blockscene.children[1].plots[4]

# %%
