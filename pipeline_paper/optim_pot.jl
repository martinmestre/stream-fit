"""Perform optimization for any fermion mass (ϵ)."""

using Pkg
Pkg.activate(".")
using PyCall
using FiniteDiff
using DelimitedFiles
using CSV
using DataFrames, DataFramesMeta
# includet("MyGIL.jl")
# using .MyGIL
using NOMAD
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

function gen_closure(g, p)
   f(x) = g(x, p)
   return f
end

function bb(x, p)
    χ² = χ²Full(x, p)
    c = -1.0
    success = true
    count_eval = true
    bb_outputs = [χ²;c]
    return (success, count_eval, bb_outputs)
end

"""Worker function."""
function worker(m, ic, r☼, lb, ub)
    local p
    len = length(lb)
    x₀ = 0.5*ones(len)
    p = (m, ic, r☼, lb, ub)
    bb_c(x) = bb(x,p)
    prob = NomadProblem(3, 2, ["OBJ", "EB"], bb_c, lower_bound=zeros(len), upper_bound=ones(len))
    result = solve(prob, x₀)
    return result
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
function cooperative(m, ic, r☼, lb_g, ub_g, n_grid)
    lb_a, ub_a, x₀_a = build_grid(lb_g, ub_g, n_grid)
    n_full = n_grid^3
    res = Vector{Float64}(undef, n_full)
    Threads.@threads for i in eachindex(x₀_a)
        println("i=$i $(Threads.threadid())")
        println("lb -- ub = , $(lb_a[i]) -- $(ub_a[i])")
        res[i] = worker(m, ic, r☼, lb_a[i], ub_a[i])
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
const maxiters = 500
const n_grid = 2
@show m sol_file r☼ maxiters

"""Running."""
# lb_l = [40.3, 29.3, 0.001]
# ub_l = [40.7, 29.7, 0.0016]
# res = worker(m, ic, r☼, lb_l, ub_l)
res = cooperative(m, ic, r☼, lb_g, ub_g, n_grid)

@show res
@show typeof(res)
# writedlm(sol_file, back_orig(sol.u, lb_l, ub_l))


