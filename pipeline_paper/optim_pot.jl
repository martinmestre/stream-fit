"""Perform optimization for ϵ ∈ [56, 370] sequentially."""

using Pkg
Pkg.activate(".")
using PyCall
using Optimization, OptimizationNOMAD
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
function χ²Full(x, p)
    θ = x[1]
    ω = x[2]
    β = x[3]/10^5
    m = p[1]
    ic = p[2]
    r☼ = p[3]
    return stream.chi2_full(θ, ω, β, m, ic, r☼)
end


"""Main function."""
function main(m, ic, r☼)
    lb = [35., 26., 1.0]
    ub = [37., 28., 1.5]
    x₀ = [36., 27., 1.25]
    p = [m, ic, r☼]
    # sol = Evolutionary.optimize(χ²Full, x₀, DE(populationSize=100),
    #                             Evolutionary.Options(iterations=100,
    #                                                  abstol=5.e-5,
    #                                                  reltol=5.e-5,
    #                                                  parallelization=:thread))
    prob = OptimizationProblem(χ²Full, x₀, p, lb=lb, ub=ub)
    sol = Optimization.solve(prob, NOMADOpt(); reltol=5.e-5)
    # open(sol_file,"a") do io
    #     println(io, "xmin = $(xmin)")
    # end
    return sol
end


# %%

"""Running."""
m = 56.0
sol = main(m, ic, r☼)
display(sol)
# %%



# %%
