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
galcen_distance = r☼*u.kpc
v_sun = astrocoords.CartesianDifferential([11.1, v_circ_sun+12.24, 7.25]*u.km/u.s)
z_sun = 0.0*u.kpc
stream_cart = astrocoords.Galactocentric(x=x*u.kpc, y=y*u.kpc, z=z*u.kpc,
                                       v_x=v_x*u.km/u.s, v_y=v_y*u.km/u.s, v_z=v_z*u.km/u.s,
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

# Accelerations
bp = galpot.PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.0,normalize=0.05)
mp = galpot.MiyamotoNagaiPotential(a=3.0/8.0,b=0.28/8.0,normalize=0.6)
hp = galpot.TriaxialNFWPotential(a=16.0/8.0,b=1.0,c=0.82,normalize=0.59)

mw = bp+mp+hp
v = Vector{Float64}(undef,3)
a_R = Vector{Float64}(undef,3)
a_z = Vector{Float64}(undef,3)
vo = 220.0u"km/s"
ro = 8.0u"kpc"
for i in 1:3
    # mw[i].turn_physical_on()
    v[i] = mw[i].vcirc(r☼/8.0)*ustrip(vo)
    a_R[i] = mw[i].Rforce(14.0/8.0, 0.0)
    a_z[i] = mw[i].zforce(14.0/8.0, 0.0)
    @show v[i] a_R[i]*galconv.force_in_kmsMyr(220.,8.) a_z[i]*galconv.force_in_kmsMyr(220.,8.)
end

v_c = sqrt(v'v)
a = [sqrt(a_R'a_R), sqrt(a_z'a_z)].*galconv.force_in_kmsMyr(220.,8.)
@show a sqrt(a'a) v_c
a = [sqrt(a_R'a_R), sqrt(a_z'a_z)]*vo^2/ro
a = uconvert.(u"km/s/Myr", a)
@show a sqrt(a'a) v_c


@show a_rar = stream.accel_mw(pot_list, stream_cart.x[50], stream_cart.y[50], stream_cart.z[50])u"(km/s)^2/kpc"
uconvert.(u"km/s/Myr",a_rar)
