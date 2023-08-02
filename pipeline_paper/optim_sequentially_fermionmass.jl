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
# using AlgebraOfGraphics, CairoMakie  # only use when needed.
using PyCall
using Optimization, OptimizationNLopt
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
# %%

"""Initial parameters, initial orbit conditions and output file."""

const param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
const θ₀, ω₀, β₀ = readdlm(param_file)
const ϵ₀ = 56.0
const E = range(ϵ₀, 370.0, step=0.1)

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)
const r☼ = 8.122

const sol_file = "param_optim_sequentially_fermionmass.txt"
const polished_file = "param_optim_sequentially_fermionmass_polished.txt"

# %%

"""Loop in ϵ."""

"""χ² wrap."""
function χ²Full(x, p)
    θ = x[1]
    ω = x[2]
    β = x[3] * 1.e-5
    ϵ = p[1]
    return stream.chi2_full(θ, ω, β, ϵ, ic, r☼)
end

function χ²Stream(x, p)
    θ = x[1]
    ω = x[2]
    β = p[1] * 1.e-5
    ϵ = p[2]
    return stream.chi2_full(θ, ω, β, ϵ, ic, r☼)
end

function χ²Core(x, p)
    θ = p[1]
    ω = p[2]
    β = x[1] * 1.e-5
    ϵ = p[3]
    return stream.chi2_full(θ, ω, β, ϵ, ic, r☼)
end

"""Optimizing sequentially by increasing the fermion mass smoothly and using
as initial parameter condition the solution of the previous mass case."""
function main(E, θ₀, ω₀, β₀, ic, r☼)
    θ, ω, βᵣ = θ₀, ω₀, β₀*10^5
    for ϵ ∈ E
        open(sol_file, "a") do f
            p = [ϵ]
            x₀ = [θ, ω, βᵣ]
            lb = [0.995θ, 0.995ω, 0.95βᵣ]
            ub = [1.005θ, 1.005ω, 1.05βᵣ]
            @show(x₀)
            prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
            sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
            x₀ = sol.u
            χ² = sol.minimum
            θ = x₀[1]
            ω = x₀[2]
            βᵣ = x₀[3]
            χ²s = stream.chi2_stream(θ, ω, βᵣ/10^5, ϵ, ic, r☼)
            χ²c = stream.chi2_core(θ, ω, βᵣ/10^5, ϵ)
            println(f, ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
            println(ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
        end
    end
end
# %%
"""Very long run."""
# main(E, θ₀, ω₀, β₀, ic, r☼)
# %%

"""Polishing for a handfull of fermion mass."""
function polish(E, in_file, out_file)
    mat = readdlm(in_file)
    df=DataFrame(ϵ=mat[:,1],θ=mat[:,2],ω=mat[:,3],βᵣ=mat[:,4]*10^5,χ²=mat[:,5],χ²s=mat[:,6],χ²c=mat[:,7])
    for ϵ ∈ E
        open(out_file, "a") do f
            df_s = @subset(df, :ϵ .== ϵ)
            θ, ω, βᵣ = df_s.θ[1], df_s.ω[1], df_s.βᵣ[1]
            p = [ϵ]
            x₀ = [θ, ω, βᵣ]
            @show x₀
            lb = [0.995θ, 0.995ω, 0.9βᵣ]
            ub = [1.005θ, 1.005ω, 1.1βᵣ]
            prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
            sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-6)
            x₀ = sol.u
            χ² = sol.minimum
            θ = x₀[1]
            ω = x₀[2]
            βᵣ = x₀[3]
            χ²s = stream.chi2_stream(θ, ω, βᵣ/10^5, ϵ, ic, r☼)
            χ²c = stream.chi2_core(θ, ω, βᵣ/10^5, ϵ)
            println(f, ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
            println(ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
        end
    end
end

"""Polishing for a handfull of fermion mass values using a two steps approach."""
function polish_in2steps(E, in_file, out_file)
    mat = readdlm(in_file)
    df=DataFrame(ϵ=mat[:,1],θ=mat[:,2],ω=mat[:,3],βᵣ=mat[:,4]*10^5,χ²=mat[:,5],χ²s=mat[:,6],χ²c=mat[:,7])
    for ϵ ∈ E
        open(out_file, "a") do f
            df_s = @subset(df, :ϵ .== ϵ)
            θ, ω, βᵣ = df_s.θ[1], df_s.ω[1], df_s.βᵣ[1]
            p = [θ, ω, ϵ]
            x₀ = [βᵣ]
            @show x₀
            lb = [0.9βᵣ]
            ub = [1.1βᵣ]
            prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
            sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
            x₀ = sol.u
            χ² = sol.minimum
            βᵣ = x₀[1]
            χ²s = stream.chi2_stream(θ, ω, βᵣ/10^5, ϵ, ic, r☼)
            χ²c = stream.chi2_core(θ, ω, βᵣ/10^5, ϵ)
            println(f, "step1", ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
            println("step1", ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)

            p = [βᵣ, ϵ]
            x₀ = [θ, ω]
            @show x₀
            lb = [0.99θ, 0.99ω]
            ub = [1.01θ, 1.01ω]
            prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
            sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
            x₀ = sol.u
            χ² = sol.minimum
            θ = x₀[1]
            ω = x₀[2]
            χ²s = stream.chi2_stream(θ, ω, βᵣ/10^5, ϵ, ic, r☼)
            χ²c = stream.chi2_core(θ, ω, βᵣ/10^5, ϵ)
            println(f,"step2", ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)
            println("step2", ϵ, "  ", θ, "  ", ω, "  ", βᵣ/10^5, "  ", χ², "  ",  χ²s, "  ", χ²c)

        end
    end
end

# %%

Eₚ = [60.0]
# polish_in2steps(Eₚ, sol_file, polished_file)

# %%
