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
"""Perform optimization for ϵ ∈ [56, 370] sequentially."""

using Pkg
Pkg.activate(".")
using PyCall
using GCMAES
using MPI
using DelimitedFiles
using CSV
using DataFrames, DataFramesMeta
# %%

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
u = pyimport("astropy.units")
importLib.reload(stream)
importLib.reload(potentials)
# %%

"""Initial orbit conditions and output file."""

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)
const r☼ = 8.122

const sol_file = "param_optim_pot.txt"


# %%

"""Loop in ϵ."""

"""χ² wrap."""
function χ²Full(x)
    θ = x[1]
    ω = x[2]
    β = x[3]/10^5
    return stream.chi2_full(θ, ω, β, m, ic, r☼)
end


"""Main function."""
function main()
    lb = [36, 25, 1.0]
    ub = [40, 30, 1.5]
    x₀ = (lb+ub)/2
    σ₀ = 0.1
    xmin, fmin, status = @mpirun GCMAES.minimize(χ²Full, x₀, σ₀, lb, ub, maxiter = 300)
    open(sol_file,"a") do io
        println(io, "xmin = $(xmin)")
    end
    println("xmin = $(xmin) with fmin = $(fmin)")
    println("status=$status")
end


# %%

"""Running."""
m = 56.0
main()
# %%



# %%
