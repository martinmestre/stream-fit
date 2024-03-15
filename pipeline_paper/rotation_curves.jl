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
"""Rotation curves computation and plotting."""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using Colors
using PyCall
using DelimitedFiles
using CSV
using DataFrames
using GLM
include("wongcolors.jl")
# %%

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
galpy_pots =  pyimport("galpy.potential")
u = pyimport("astropy.units")
importLib.reload(stream)
importLib.reload(potentials)

# %%
# Parameters and initial conditions.
ϵ = 56
param_file = "serafin/sol_optim_pot_m$(Int(ϵ)).txt"
θ, ω, β = vec(readdlm(param_file))
W = θ+ω
param = [ϵ, θ, W, β]


ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122


# %%
# Fermionic-MW solution and rotation curve
pot_list = stream.pot_model(ϵ, θ, W, β)
halo = pot_list[4]
r = 10 .^ range(log10(halo.r_s[begin]),log10(100.),length=2000)
v_ferm = [stream.rot_vel_mw(pot_list, x) for x in r]



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
df_nfw = DataFrame(grp="NFW-MW", r=r, v=v_nfw)
df_ferm = DataFrame(grp="Fermionic-MW", r=r, v=v_ferm)
df_model = append!(df_nfw, df_ferm)


 # Reading rotation curve observables
 # Set rotation data

df_Sof13 = DataFrame(CSV.File("data_rotcurves/vel_Sofue13.txt"; delim=" ", ignorerepeated=true))
df_Sof20 = DataFrame(CSV.File("data_rotcurves/vel_Sofue20.txt"; delim=" ", ignorerepeated=true))
df_Eilers = DataFrame(CSV.File("data_rotcurves/vel_Eilers.txt"; delim=" ", ignorerepeated=true))
df_Ou = DataFrame(CSV.File("data_rotcurves/vel_Ou23.txt"; delim=" ", ignorerepeated=true))
df_Jiao = DataFrame(CSV.File("data_rotcurves/vel_Jiao23.txt"; delim=" ", ignorerepeated=true))

df_Sof13.r = df_Sof13.r/1.0e3
df_Sof13.err_r =  df_Sof13.err_r/1.0e3
rename!(df_Sof20, [:Radius, :Velocity, :Error] .=> [:r, :v, :err_v])
insertcols!(df_Sof13, 1, :grp=>fill("Sofue 2013",nrow(df_Sof13)))
insertcols!(df_Sof20, 1, :grp=>fill("Sofue 2020",nrow(df_Sof20)))
insertcols!(df_Eilers, 1, :grp=>fill("Eilers 2019",nrow(df_Eilers)))
insertcols!(df_Ou, 1, :grp=>fill("Ou 2023",nrow(df_Ou)))
insertcols!(df_Jiao, 1, :grp=>fill("Jiao 2023",nrow(df_Jiao)))
@show names(df_Sof13) names(df_Sof20) names(df_Eilers)
@show names(df_Ou) names(df_Jiao)

df_obs = vcat(df_Sof20, df_Eilers, df_Jiao, cols=:union)
for col in eachcol(df_obs)
    replace!(col,missing=>0)
end
# %%

# Paper plot
labels = ["Fermionic-MW","NFW-MW"]
lw = 4

col = wongcolors()
pal = (color=[col[4], col[2], :black, col[3], col[5]],
    patchcolor=wongcolors(),
    marker=[:utriangle, :utriangle, :circle, :rect, :diamond, :dtriangle, :pentagon, :xcross],
    linestyle=[:solid, :dash, :dashdot, :dot, :dashdotdot],
    side=[:left, :right],
)

set_aog_theme!()
update_theme!(palette=pal,Axis=(topspinevisible=true, rightspinevisible=true,
topspinecolor=:darkgray, rightspinecolor=:darkgray,
xticksmirrored = true, yticksmirrored = true))
size_inches = (6.2*2, 4*2)
size_pt = 72 .* size_inches
fig = Figure(resolution = size_pt, fontsize = 37)
gridpos = fig[1, 1]

rv_map = mapping(:r => L"r~[\textrm{kpc}]", :v=> L"v_{\textrm{circ}}~[\textrm{km s}^{-1}]";
                 color = :grp => "", marker = :grp => "", linestyle = :grp => "")
rv_ev_map = mapping(:r=> L"r~[\textrm{kpc}]",:v=> L"v_{\textrm{circ}}~[\textrm{km s}^{-1}]",
                  :err_v; color =   :grp => "", marker = :grp => "", linestyle = :grp=>"")
rv_ev_ud_map = mapping(:r=> L"r~[\textrm{kpc}]",:v=> L"v_{\textrm{circ}}~[\textrm{km s}^{-1}]",
                    :e_down, :e_up; color =   :grp => "", marker = :grp => "", linestyle = :grp=>"")
rv_er_map = mapping(:r=> L"r~[\textrm{kpc}]",:v=> L"v_{\textrm{circ}}~[\textrm{km s}^{-1}]",
                    :err_r; color =   :grp => "", marker = :grp => "", linestyle = :grp=>"")

plt_model = data(df_model) * rv_map * visual(Lines, linewidth=lw)
plt_obs   = data(df_obs)   * rv_map * visual(Scatter)
plt_ev    = data(df_obs)   * rv_ev_map * visual(Errorbars)
plt_ev_ud =  data(df_obs)   * rv_ev_ud_map * visual(Errorbars)
plt_er = data(df_obs)   * rv_er_map * visual(Errorbars;direction=:x)


plt = plt_model+plt_obs+plt_ev+plt_ev_ud
f=draw!(gridpos, plt, axis=(;limits=((0,42),(80,270)),
    xgridvisible=false, ygridvisible=false))
ax = f[1,1].axis
vspan!(ax,11.539089193812874,16.35356513868187,color=(:black,0.15))

legend!(gridpos, f; tellwidth=false, halign=:right, valign=:bottom,
        margin=(10, 10, 10, 10), patchsize=(50,20), nbanks=2, framevisible=false, labelsize=25)

# Lines re-styling
col = wongcolors()
lineas = fig.content[1].scene.plots
println("lineas=$lineas")


# Legend re-styling
leg = fig.content[2]
_lines = leg.blockscene.children[1].plots
for l in _lines
    l.linewidth = 5
end

println("_lines=$_lines")
deleteat!(_lines,[2,4,6,7,8,10,12,13,14,16])
display(fig)
save("paper_plots/rotation_curves.pdf", fig, pt_per_unit = 1)
println("first plot done.")
# %%

# Zoomed in plot (not in paper)
fig = Figure(resolution = size_pt, fontsize = 37)
gridpos = fig[1, 1]
f=draw!(gridpos, plt, axis=(;limits=((5,20),(220,270)),
    xgridvisible=false, ygridvisible=false))
ax = f[1,1].axis
vspan!(ax,11.539089193812874,16.35356513868187,color=(:black,0.15))
legend!(gridpos, f; tellwidth=false, halign=:left, valign=:top,
        margin=(10, 10, 10, 10), patchsize=(30,20), nbanks=2, framevisible=true, labelsize=25)
leg = fig.content[2]
_lines = leg.blockscene.children[1].plots
for l in _lines
    l.linewidth = 4
end
println("_lines=$_lines")
deleteat!(_lines,[3,4,5,6,8,9,12,13,15,16])
display(fig)
save("paper_plots/rotation_curves_zoom.pdf", fig, pt_per_unit = 1)
println("Second plot done.")



# %%
# Fit data of the fermionic rotation curve.
bool = 14.5 .< df_ferm.r .< 26.5
mydata = DataFrame(x=df_ferm.r[bool],y=df_ferm.v[bool])
ols = lm(@formula(y ~ 1+x), mydata)


# %%
fig.content[2].blockscene.children[1].plots



# Tests...
# %%
n = 20
df_a = DataFrame(grp=fill("a", n), x=collect(1:n), y=fill(2, n))
df_b = DataFrame(grp=fill("b", n), x=collect(1:n), y=fill(3, n), z=fill(0, n))

set_aog_theme!()
size_inches = (6.2*2, 3*2)
size_pt = 72 .* size_inches
fig = Figure(resolution = size_pt, fontsize = 30)
gridpos = fig[1, 1]
plt_a = data(df_a) *
    mapping(:x => L"x", :y=> L"y"; color = :grp => "", marker=:grp=> "", linestyle = :grp => "") *
    visual(Lines, linewidth=lw)
plt_b = data(df_b) *(
    mapping(:x => L"x", :y=> L"y"; color = :grp => "", marker=:grp=>"", linestyle= :grp => "") *
    visual(Scatter))
    # + mapping(:x => L"x",:y=> L"y",:z;color=:grp=> "", marker=:grp=> "", linestyle = :grp => "")*visual(Errorbars))
f = draw!(gridpos, plt_a+plt_b, axis=(;limits=((0,n),(0,4)),
    xgridvisible=false, ygridvisible=false))
legend!(gridpos, f; tellwidth=false, halign=:center, valign=:bottom, margin=(10, 10, 10, 10), patchsize=(20,10))

display(fig)
leg = fig.content[2]
leg_items = leg.blockscene.children[1].plots
display(leg_items)
# %%
