"""Perform optimization for ϵ ∈ [56, 370] sequentially."""

using Pkg
Pkg.activate(".")
using PyCall
using Optimization, OptimizationOptimisers
using FiniteDiff
using DelimitedFiles
using CSV
using DataFrames, DataFramesMeta
# includet("MyGIL.jl")
# using .MyGIL

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
    @show  θ, ω, β
    return stream.chi2_full(θ, ω, β, m, ic, r☼)
end


function χ²Full_parallel(x, p)
    println("sizes=",size(x,1), size(x,2))
    len = size(x,1)
    fitness = zeros(len)
    Threads.@threads for i in 1:size(x,1)
        MyGIL.pylock() do
            @show χ²Full(x[i,:], p)
            fitness[i] = χ²Full(x[i,:], p)
        end
    end
    return fitness
end

"""Worker function."""
function worker(m, ic, r☼, maxiters, lb, ub)
    len = length(lb)
    x₀ = 0.5*ones(len)
    p = (m, ic, r☼, lb, ub)
    prob = OptimizationProblem(χ²Full_parallel, x₀, p, lb=zeros(len), ub=ones(len))
    println("prob=$prob")
    sol = Optimization.solve(prob, ECA(); maxiters=maxiters, reltol=5.e-7, parallel_evaluation=true, use_initial=true)
    return sol
end

"""Build grid."""
function build_grid(lb_g, ub_g, n_grid)
    n_full = n_grid^3
    lb_a = Vector{Vector{Float64}}(undef,n_full)
    ub_a = Vector{Vector{Float64}}(undef,n_full)
    x₀_a = Vector{Vector{Float64}}(undef,n_full)
    c₁ = collect(range(lb_g[1], ub_g[1], n_grid+1))
    c₂ = collect(range(lb_g[2], ub_g[2], n_grid+1))
    c₃ = collect(range(lb_g[3], ub_g[3], n_grid+1))
    for i ∈ 1:n_grid
        for j ∈ 1:n_grid
            for k ∈ 1:n_grid
                n = (i-1)*n_grid^2+(j-1)*n_grid+k
                lb_a[n] = [c₁[i], c₂[j], c₃[k]]
                ub_a[n] = [c₁[i+1], c₂[j+1], c₃[k+1]]
                x₀_a[n] = 0.5*(lb_a[n]+ub_a[n])
            end
        end
    end
    return lb_a, ub_a, x₀_a
end

"""Parallel function."""
function cooperative(m, ic, r☼, maxiters, lb_g, ub_g, n_grid)
    lb_a, ub_a, x₀_a = build_grid(lb_g, ub_g, n_grid)
    Threads.@threads for i in eachindex(x₀_a)
        println("i=$i $(Threads.threadid())")
        println("lb_ub = , $(lb_a[i]) -- $(ub_a[i])")
        worker(m, ic, r☼, maxiters, lb_a[i], ub_a[i])
    end
end
# %%

"""Initial orbit conditions file."""

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)

"""Metaparameters."""
const m = 300.0
const sol_file = "param_optim_pot_m$(Int(m)).txt"
const r☼ = 8.122
const lb_g = [35., 25., 1.e-5]
const ub_g = [45., 31., 0.005]
const maxiters = 1000
const n_grid = 2
@show m sol_file r☼ maxiters

"""Running."""
lb_l = [40., 29., 0.001]
ub_l = [41., 30., 0.002]
# xᵢ=[40.654391903440164, 29.52500686293749, 0.0013010438131347817]
# lb_l = xᵢ-[0.01, 0.01, 0.0001]
# ub_l = xᵢ+[0.01, 0.01, 0.0001]
# sol = worker(m, ic, r☼, maxiters, lb_l, ub_l)
len = length(lb_l)
x₀ = 0.5*ones(len)
p = (m, ic, r☼, lb_l, ub_l)
optfun = OptimizationFunction(χ²Full, AutoFiniteDiff())
prob = OptimizationProblem(χ²Full, x₀, p)
println("prob=$prob")
sol = Optimization.solve(prob, Optimisers.Descent())
@show sol
writedlm(sol_file, back_orig(sol.u, lb_l, ub_l))
# %%
# cooperative(m, ic, r☼, maxiters, lb_g, ub_g, n_grid)



