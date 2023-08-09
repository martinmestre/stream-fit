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
using Evolutionary
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
println("threads=", Threads.nthreads())

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
    lb = [35, 37, 1.0]
    ub = [25, 30, 1.5]
    x₀ = [36., 27., 1.25]
    res = Evolutionary.optimize(χ²Full, x₀, DE(populationSize=100),
                                Evolutionary.Options(iterations=100,
                                                     abstol=5.e-5,
                                                     reltol=5.e-5,
                                                     parallelization=:thread))
    # open(sol_file,"a") do io
    #     println(io, "xmin = $(xmin)")
    # end
    return res
end


# %%

"""Running."""
m = 56.0
res = main()
display(res)
# %%



# %%
