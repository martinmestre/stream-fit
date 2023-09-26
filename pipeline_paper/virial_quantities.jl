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

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)

# Parameters and initial conditions.
param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
θ, Δθ, β = readdlm(param_file)
ϵ = 56.0
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
