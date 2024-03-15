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
"""Compute virial quantities of our solution."""

using Pkg
Pkg.activate(".")
using PyCall
using CSV
using DelimitedFiles
using Unitful, UnitfulAstro

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)

# Parameters and initial conditions.
ϵ = 56.0
param_file = "serafin/sol_optim_pot_m$(Int(ϵ)).txt"
θ, Δθ, β = vec(readdlm(param_file))
W = θ+Δθ
param = [ϵ, θ, W, β]

ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122
# %%

pot_list = stream.pot_model(ϵ, θ, W, β)
rfk = pot_list[4]
@show rfk.r_max==rfk.r_s[end]
r_vir = rfk.r_max
m_vir = rfk.mass_wrap(r_vir)[1]
@show r_vir m_vir;
# %%

# Data from table 3 in Gibbons et al. 2014
mass_mw = 2.9e11
σ = 0.4e11
σx2 = 0.9e11
unit_mass_Allen = 2.32e7
mass_disks = 2.0*1700.0*unit_mass_Allen
mass_bulge = 460.0*unit_mass_Allen
mass_barions = mass_disks+mass_bulge
mass_halo = mass_mw - σx2 - mass_barions
println("Disks mass = ", mass_disks/1.e11, " x 10^11 M_⊙")
println("Barions mass = ", mass_barions/1.e11, " x 10^11 M_⊙")
println("Halo mas = ", mass_halo/1.e11, " x 10^11 M_⊙")
@show m_vir+mass_barions

# %%
# Local density
ρ☼ = rfk.rho_spl(r☼)[1]
@show ρ☼
# fact value is the results uconvert(u"GeV/cm^3",1.0unit(ρ☼)*u"c"^2)
fact = 3.7965164572915704e-8
@show fact*ρ☼


# %%
# Enclosed mass
r☼ = 8.122u"kpc"
r_p_as = 0.01427*u"arcsecond"  # New Edward values
r_a_as = 0.23623*u"arcsecond"
r_p = uconvert(u"kpc", r_p_as*r☼)
r_a = uconvert(u"kpc", r_a_as*r☼)
@show r_p r_a
m_p = rfk.mass_wrap(r_p/u"kpc")[1]
Δm = rfk.mass_wrap(r_a/u"kpc")[1]-rfk.mass_wrap(r_p/u"kpc")[1]
m_core = 3.5e6
@show m_p Δm/m_core

m12 = rfk.mass_wrap(12.0)[1]
m40 = rfk.mass_wrap(40.0)[1]
m80 = rfk.mass_wrap(80.0)[1]

@show m12 m40 m80;

# %%
