#!/usr/bin/env julia

"""Perform optimization for any fermion mass (ϵ) by running in parallel
for a grid in a global region of paramter space.
Using Distributed.jl
"""

using Pkg
Pkg.activate(".")
using Distributed
using SlurmClusterManager
addprocs(SlurmManager(;launch_timeout=300.0))


@show nprocs()
@show nworkers()

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

@everywhere begin
    using PythonCall
    using Optimization, OptimizationNOMAD
    using FiniteDiff
    using DelimitedFiles
    using CSV
    using DataFrames, DataFramesMeta
end

@everywhere begin
    pyimport("sys")."path".append("")
    stream = pyimport("stream")
    potentials = pyimport("potential_classes")
    u = pyimport("astropy.units")
end
    # %%
@everywhere begin
    """χ² wrap."""
    function χ²Full(x, p)
        m = p[1]
        ic = p[2]
        r☼ = p[3]
        θ, ω, β = x
        return pyconvert(Float64,stream.chi2_full(θ, ω, β, m, ic, r☼))
    end


    """Worker function."""
    function worker(i, sol_dir, m, ic, r☼, lb, ub)
        println("Now at id:$(myid()), host:$(gethostname())")
        x₀ = 0.5*(lb+ub)
        p = (m, ic, r☼)
        prob = OptimizationProblem(χ²Full, x₀, p, lb=lb, ub=ub)
        sol = Optimization.solve(prob, NOMADOpt(); display_degree=0, maxiters=700)
        worker_file = "$(sol_dir)/worker_optim_pot_m$(Int(m))_i$i.txt"
        worker_sol = ("Minimizer = $(sol.u)", "Minimum = $(sol.objective)")
        writedlm(worker_file, worker_sol)
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
    function cooperative(sol_dir, m, ic, r☼, lb_g, ub_g, n_grid)
        lb_a, ub_a, x₀_a = build_grid(lb_g, ub_g, n_grid)
        pars(i) = [i, sol_dir, m, ic, r☼, lb_a[i], ub_a[i]]
        res = pmap(i->worker(pars(i)...), eachindex(x₀_a))
        return res
    end
    # %%
end


"""Initial orbit conditions file."""

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = vec(readdlm(ic_file))

"""Metaparameters."""
const i_m = parse(Int,ARGS[1])
const m_a = [56., 100., 200., 300., 360., 378.]
const m = m_a[i_m]
const sol_dir = "sol_dir_optim_pot_m$(Int(m))"
const sol_file = "sol_optim_pot_m$(Int(m)).txt"
const r☼ = 8.122

const lb_g = [[35., 26., 1.0e-5], [36., 27., 1.2e-5], [37., 28., 5.0e-5],
              [38., 29., 3.5e-4], [40., 29., 1.3e-3], [43., 29.6, 3.0e-3]]
const ub_g = [[39., 30., 1.5e-5], [40., 31., 1.0e-4], [41., 32., 1.0e-3],
              [42., 32., 3.0e-3], [44., 32., 4.0e-3], [47., 36., 1.0e-2]]
const n_grid = 20
@show m sol_file r☼ lb_g ub_g

if !isdir(sol_dir)
    run(`mkdir $sol_dir`)
end

# """Running."""
sol = cooperative(sol_dir, m, ic, r☼, lb_g[i_m], ub_g[i_m], n_grid)
@show sol
obj = [sol[i].objective for i ∈ eachindex(sol)]
min, index = findmin(obj)
best_u = sol[index].u
best = ("Minimizer = $(best_u)", "Minimum = $(min)")
writedlm(sol_file, best)
# %%




