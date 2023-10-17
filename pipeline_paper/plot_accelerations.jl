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

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)
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
ϕ₂ = fϕ₂(ϕ₁ₒ)
d☼ = fd☼(ϕ₁ₒ)
μ_ra = fμ_ra(ϕ₁ₒ)
μ_dec = fμ_dec(ϕ₁ₒ)
v☼ = fv☼(ϕ₁ₒ)

# Malhan (MWPotential2014) solution.
temp = readdlm(orbit_nfw_file)
boolϕ₁ = (temp[1,:].>-95 .&& temp[1,:].<15)
fϕ₂, fd☼, fμ_ra, fμ_dec, fv☼ = [CubicSpline(temp[1,:][boolϕ₁], temp[i,:][boolϕ₁]) for i=2:6]
ϕ₂ₘ = fϕ₂(ϕ₁ₒ)
d☼ₘ = fd☼(ϕ₁ₒ)
μ_raₘ = fμ_ra(ϕ₁ₒ)
μ_decₘ = fμ_dec(ϕ₁ₒ)
v☼ₘ = fv☼(ϕ₁ₒ);
# %%

