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
using LaTeXStrings
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
temp = stream.orbit_model_ext(ic..., pot_list, r☼)
# Cubic spline for solution:
boolϕ₁ = (temp[1].>-95 .&& temp[1].<15)
fϕ₂, fd☼, fμ_ra, fμ_dec, fv☼ = [CubicSpline(temp[1][boolϕ₁], temp[i][boolϕ₁]) for i=2:6]
fx = CubicSpline(temp[1][boolϕ₁], temp[7][boolϕ₁])
fy = CubicSpline(temp[1][boolϕ₁], temp[8][boolϕ₁])
fz = CubicSpline(temp[1][boolϕ₁], temp[9][boolϕ₁])
fvx = CubicSpline(temp[1][boolϕ₁], temp[11][boolϕ₁])
fvy = CubicSpline(temp[1][boolϕ₁], temp[12][boolϕ₁])
fvz = CubicSpline(temp[1][boolϕ₁], temp[13][boolϕ₁])
ϕ₂ = fϕ₂(ϕ₁ₒ)
d☼ = fd☼(ϕ₁ₒ)
μ_ra = fμ_ra(ϕ₁ₒ)
μ_dec = fμ_dec(ϕ₁ₒ)
v☼ = fv☼(ϕ₁ₒ)
x, y, z = fx(ϕ₁ₒ), fy(ϕ₁ₒ), fz(ϕ₁ₒ)
v_x, v_y, v_z = fvx(ϕ₁ₒ), fvy(ϕ₁ₒ), fvz(ϕ₁ₒ)
v_circ_sun = stream.rot_vel_mw(pot_list, r☼)
galcen_distance = r☼*au.kpc
v_sun = astrocoords.CartesianDifferential([11.1, v_circ_sun+12.24, 7.25]*au.km/au.s)
z_sun = 0.0*au.kpc
stream_cart = astrocoords.Galactocentric(x=x*au.kpc, y=y*au.kpc, z=z*au.kpc,
                                       v_x=v_x*au.km/au.s, v_y=v_y*au.km/au.s, v_z=v_z*au.km/au.s,
                                       galcen_distance=galcen_distance, galcen_v_sun=v_sun, z_sun=z_sun)
# Malhan (MWPotential2014) solution.
temp = readdlm(orbit_nfw_file)
boolϕ₁ = (temp[1,:].>-95 .&& temp[1,:].<15)
fϕ₂, fd☼, fμ_ra, fμ_dec, fv☼ = [CubicSpline(temp[1,:][boolϕ₁], temp[i,:][boolϕ₁]) for i=2:6]
ϕ₂ₘ = fϕ₂(ϕ₁ₒ)
d☼ₘ = fd☼(ϕ₁ₒ)
μ_raₘ = fμ_ra(ϕ₁ₒ)
μ_decₘ = fμ_dec(ϕ₁ₒ)
v☼ₘ = fv☼(ϕ₁ₒ)
#%%

# Acceleration
function acceleration_nfw_mw(pot::PyObject, x::Vector{T})  where {T<:Number}
    a_R = Vector{Float64}(undef,3)
    a_z = Vector{Float64}(undef,3)
    R = sqrt(x[1:2]'x[1:2])
    z = x[3]
    r = [R, z]
    for i in 1:3
        a_R[i] = pot[i].Rforce(r...)[1]
        a_z[i] = pot[i].zforce(r...)[1]
    end
    return [sum(a_R), sum(a_z)]
end
function acceleration_rfk_mw(pot_list::Vector{PyObject}, x::Vector{T})  where {T<:Number}
    a = stream.accel_mw(pot_list, x...)u"(km/s)^2/kpc"
    a = uconvert.(u"km/s/Myr", a)
    eᵣ = x[1:2]/sqrt(x[1:2]'x[1:2])
    a_R = a[1:2]'eᵣ
    return [a_R, a[3]]
end

bp = galpot.PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.0,normalize=0.05)
mp = galpot.MiyamotoNagaiPotential(a=3.0/8.0,b=0.28/8.0,normalize=0.6)
hp = galpot.TriaxialNFWPotential(a=16.0/8.0,b=1.0,c=0.82,normalize=0.59)
mw = bp+mp+hp
for i in 1:3
    mw[i].turn_physical_on()
end
jᵣ = 80
x = [stream_cart.x[jᵣ], stream_cart.y[jᵣ], stream_cart.z[jᵣ]]
a_nfw = acceleration_nfw_mw(mw, x/8.0)
a_rfk = acceleration_rfk_mw(pot_list, x)
@show a_nfw./a_rfk sqrt(a_nfw'a_nfw/a_rfk'a_rfk)
@show a_nfw a_rfk;
x = [[stream_cart.x[i], stream_cart.y[i], stream_cart.z[i]] for i in eachindex(stream_cart.x)]
R = [sqrt(stream_cart.x[i]^2+stream_cart.y[i]^2) for i in eachindex(stream_cart.x)]
z = [stream_cart.z[i] for i in eachindex(stream_cart.x)]
a_nfw=[acceleration_nfw_mw(mw, x[i]/8.0) for i in eachindex(x)]
a_rfk=[acceleration_rfk_mw(pot_list, x[i]) for i in eachindex(x)]
#%%

df = DataFrame(ϕ₁=ϕ₁ₒ,
                R=R,
                z=z,
                a_nfw_R=[a_nfw[i][1] for i ∈ eachindex(a_nfw)],
                a_nfw_z=[a_nfw[i][2] for i ∈ eachindex(a_nfw)],
                a_rfk_R=ustrip.([a_rfk[i][1] for i ∈ eachindex(a_rfk)]),
                a_rfk_z=ustrip.([a_rfk[i][2] for i ∈ eachindex(a_rfk)])
)
#%%
set_aog_theme!()
update_theme!(Axis=(topspinevisible=true, rightspinevisible=true,
topspinecolor=:darkgray, rightspinecolor=:darkgray,
xticksmirrored = true, yticksmirrored = true))
size_inches = (6.2*2, 3*2)
size_pt = 72 .* size_inches
lw = 5
fig = Figure(resolution = size_pt, fontsize = 33)
gridpos = fig[1, 1]
labels = [L"radial (NFW-MW)", L"z (NFW-MW)", L"radial (Fermionic-MW)", L"z (Fermionic-MW)"]
grp = dims(1) => renamer(labels) => ""
plt = data(df) *
    mapping(:ϕ₁ => L"ϕ_1~[°]", [4,5,6,7] .=> L"a~[\mathrm{km~s}^{-1}~\mathrm{Myr}^{-1}]";
        color = grp,
        linestyle = grp
    ) *
    visual(Lines, linewidth=lw)
f = draw!(gridpos, plt, axis=(;limits=((-90,10),(-4, 1)),
    xgridvisible=false, ygridvisible=false))

legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(55,40))
display(fig)
#%%

set_aog_theme!()
update_theme!(Axis=(topspinevisible=true, rightspinevisible=true,
topspinecolor=:darkgray, rightspinecolor=:darkgray,
xticksmirrored = true, yticksmirrored = true))
size_inches = (6.2*2, 3*2)
size_pt = 72 .* size_inches
lw = 5
fig = Figure(resolution = size_pt, fontsize = 33)
gridpos = fig[1, 1]
labels = ["R (NFW-MW)", "z (NFW-MW)", "R (Fermionic-MW)", "z (Fermionic-MW)"]
grp = dims(1) => renamer(labels) => ""
plt = data(df) *
    mapping(:R => L"R~[kpc]", [4,5,6,7] .=> L"a~[\mathrm{km s}^{-1} \mathrm{Myr}^{-1}]";
        color = grp,
        linestyle = grp
    ) *
    visual(Lines, linewidth=lw)
f = draw!(gridpos, plt, axis=(xgridvisible=false, ygridvisible=false))

legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(55,40))
display(fig)
#%%

set_aog_theme!()
update_theme!(Axis=(topspinevisible=true, rightspinevisible=true,
topspinecolor=:darkgray, rightspinecolor=:darkgray,
xticksmirrored = true, yticksmirrored = true))
size_inches = (6.2*2, 3*2)
size_pt = 72 .* size_inches
lw = 5
fig = Figure(resolution = size_pt, fontsize = 33)
gridpos = fig[1, 1]
labels = ["R (NFW-MW)", "z (NFW-MW)", "R (Fermionic-MW)", "z (Fermionic-MW)"]
grp = dims(1) => renamer(labels) => ""
plt = data(df) *
    mapping(:z => L"z~[kpc]", [4,5,6,7] .=> L"a~[\mathrm{km s}^{-1} \mathrm{Myr}^{-1}]";
        color = grp,
        linestyle = grp
    ) *
    visual(Lines, linewidth=lw)
f = draw!(gridpos, plt, axis=(xgridvisible=false, ygridvisible=false))

legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(55,40))
display(fig)
#%%

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
    visual(Lines, linewidth=lw)
f = draw!(gridpos, plt, axis=(xgridvisible=false, ygridvisible=false, limits=((0,16),(0, 16))))
legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(55,40))
us = df.a_rfk_R
vs = df.a_rfk_z
strength = @. sqrt(us^2+vs^2)
@. us = us/strength
@. vs = vs/strength
arrows!(df.R, df.z, df.a_rfk_R, df.a_rfk_z, arrowsize = strength, lengthscale = 1)
display(fig)


