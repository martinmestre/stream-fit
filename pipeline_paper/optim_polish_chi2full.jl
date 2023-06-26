"""Optimization trials."""

# %%
using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using PyCall
using Optimization
# using OptimizationNLopt
# using OptimizationCMAEvolutionStrategy
# OES = OptimizationCMAEvolutionStrategy
using Evolutionary
using DelimitedFiles


# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")

# %%
function χ²Core(x, p)
    β = x[1] * 1.e-5
    θ = p[1]
    Δθ = p[2]
    ϵ = p[3]
    return stream.chi2_core(θ, Δθ, β, ϵ)
end

function χ²Stream(x, p)
    θ = x[1]
    Δθ = x[2]
    β = p[1] * 1.e-5
    ϵ = p[2]
    println("inside chi2stream x=$x")
    return stream.chi2_stream(θ, Δθ, β, ϵ, ic, r☼)
end

function χ²Full(x)
    θ = x[1]
    Δθ = x[2]
    β = x[3] * 1.e-5
    # ϵ = p[1]
    println("inside chi2full x=$x")
    return stream.chi2_full(θ, Δθ, β, ϵ, ic, r☼)
end

function main()
    # %%
    infile = "param_fit_pot-slice_from_IbataPolysGaiaDR2-data.txt"
    outfile = "param_optim_polish_chi2full.txt"

    # Initial conditions from Galpy's MWPotential2014 model (param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt):
    ic = [149.16883497, 36.50694499, 7.94771853, -6.94192661, -12.48844894, -18.32838854]
    r☼ = 8.122
    # θ₀, Δθ₀ = open(readdlm, infile)
    θ₀, Δθ₀ = 37.5, 27.5
    βᵣ₀ =  1.2  # it goes without the factor 10^(-5).
    ϵ = 56.0
    # %%
    p = [ϵ]
    x₀ = [θ₀, Δθ₀, βᵣ₀]
    println("x₀=$x₀")
    # lb = [0.99θ₀, 0.99Δθ₀, 0.99βᵣ₀]
    # ub = [1.01θ₀, 1.01Δθ₀, 1.1βᵣ₀]
    lb = [35., 25., 0.9βᵣ₀]
    ub = [40., 30., 1.1βᵣ₀]
    prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
    # sol = solve(prob, CMAEvolutionStrategyOpt(), reltol=5.0e-5
    sol = Evolutionary.optimize(χ²Full,BoxConstraints(lb, ub), x₀, GA(selection=uniformranking(5),
    mutation=flip, crossover=SPX), Evolutionary.Options(iterations=10, parallelization=:thread))
    x₀ = sol.u
    χ² = sol.minimum
    println("after fit: x₀=", x₀, "χ²Stream=$χ²Stream","χ²Full=$χ²Full")
    # θ₀ = x₀[1]
    # Δθ₀ = x₀[2]
    # βᵣ₀ = x₀[3]
    # χ²s = stream.chi2_stream(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ, ic, r☼)
    # χ²c = stream.chi2_core(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ)
    # write(outfile, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5, "  ", χ², "  ",  χ²s, "  ", χ²c)
    # println("χ²:F-S-C:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5,
    #         "  ", χ², "  ",  χ²s, "  ", χ²c)
end


