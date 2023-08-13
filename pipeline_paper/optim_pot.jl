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

println("Threads=", Threads.nthreads())

"""Loop in ϵ."""

"""Anti-normalization function."""
function back_orig(x, a, b)
    return (b-a).*x + a
end

"""χ² wrap."""
function χ²Full(x, p)
    m = p[1]
    ic = p[2]
    r☼ = p[3]
    lb = p[4]
    ub = p[5]
    θ, ω, β = back_orig(x, lb, ub)
    return stream.chi2_full(θ, ω, β, m, ic, r☼)
end


"""Worker function."""
function worker(m, ic, r☼, reltol, maxiters, lb, ub)
    len = length(lb)
    x₀ = 0.5*ones(len)
    p = (m, ic, r☼, lb, ub)
    prob = OptimizationProblem(χ²Full, x₀, p, lb=zeros(len), ub=ones(len))
    @show prob
    sol = Optimization.solve(prob, NOMADOpt(); reltol=reltol, maxiters=maxiters)
    return sol
end

"""Build grid."""
function build_grid(g_lb, g_ub, n_grid)
    n_full = n_grid^3
    lb_a = Vector{Vector{Float64}}(undef,n_full)
    ub_a = Vector{Vector{Float64}}(undef,n_full)
    x₀_a = Vector{Vector{Float64}}(undef,n_full)
    c₁ = collect(range(g_lb[1], g_ub[1], n_grid))
    c₂ = collect(range(g_lb[2], g_ub[3], n_grid))
    c₃ = collect(range(g_lb[3], g_ub[3], n_grid))
        for i ∈ 1:n_grid
            for j ∈ 1:n_grid
                for k ∈ 1:n_grid
                    n = (i-1)*n_grid^2+(j-1)*n_grid+k
                    lb_a[n] = [c₁[i], c₂[j], c₃[k]]
                    x₀_a[n] = []
    return lb_a, ub_a, x₀_a
end

"""Parallel function."""
function cooperative(m, ic, r☼, reltol, maxiters, g_lb, g_ub, n_grid)
    t_lb, t_ub = build_grid(g_lb, g_ub, ngrid)
end
# %%

"""Initial orbit conditions file."""

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)

"""Metaparameters."""
const m = 56.0
const sol_file = "param_optim_pot_m$(Int(m)).txt"
const r☼ = 8.122
const glob_lb = [35., 25., 1.e-5]
const glob_ub = [45., 31., 0.005]
const reltol = 5.0e-5
const maxiters = 2
@show m sol_file r☼ reltol maxiters

"""Running."""
sol = worker(m, ic, r☼, lb, ub, reltol, maxiters)
# @show sol
# writedlm(sol_file, sol.u)
# %%



