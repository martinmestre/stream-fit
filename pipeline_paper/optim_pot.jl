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

println("threads=", Threads.nthreads())

"""Loop in ϵ."""

"""χ² wrap."""
function χ²Full(x, p)
    θ = x[1]
    ω = x[2]
    β = x[3]
    m = p[1]
    ic = p[2]
    r☼ = p[3]
    return stream.chi2_full(θ, ω, β, m, ic, r☼)
end


"""Main function."""
function main(m, ic, r☼)
    # lb = [41., 28., 0.001] # for m = 360
    # ub = [45., 31., 0.005]
    # lb = [39., 28., 0.0005]  # for m = 300
    # ub = [42., 31., 0.003 ]
    # lb = [38., 28., 0.0001] # for m = 200
    # ub = [40., 30., 0.0006]
    # lb = [36., 27., 4.0e-5] # for m=100
    # ub = [38., 29., 7.0e-5]
    lb = [34., 26., 1.0e-5] # for m=56
    ub = [38., 30., 1.5e-5]
    x₀ = 0.5*(lb+ub)
    p = [m, ic, r☼]
    prob = OptimizationProblem(χ²Full, x₀, p, lb=lb, ub=ub)
    sol = Optimization.solve(prob, NOMADOpt(); reltol=5.e-5, maxiters=1000)
    return sol
end
# %%

"""Initial orbit conditions file."""

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)
const r☼ = 8.122
# %%

"""Running."""
const m = 56.0
const sol_file = "param_optim_pot_m$(Int(m)).txt"
@show m sol_file
sol = main(m, ic, r☼)
@show sol
writedlm(sol_file, sol.u)
# %%



