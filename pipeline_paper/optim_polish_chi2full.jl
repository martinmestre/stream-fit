"""Optimization trials."""

# %%
using AlgebraOfGraphics, CairoMakie
using PyCall
using GalacticOptim
using GalacticOptimJL
using GalacticNLopt

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
    return stream.chi2_stream(θ, Δθ, β, ϵ, ic, r☼)
end

function χ²Full(x, p)
    θ = x[1]
    Δθ = x[2]
    β = x[3] * 1.e-5
    ϵ = p[1]
    return stream.chi2_full(θ, Δθ, β, ϵ, ic, r☼)
end

# %%
# Initial conditions from Galpy's MWPotential2014 model (param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt):
ic = [149.16883497, 36.50694499, 7.94771853, -6.94192661, -12.48844894, -18.32838854]
r☼ = 8.122

θ₀, Δθ₀, βᵣ₀ =  3.612769971239480782e+01, 2.739433669170758989e+01, 1.2e-5
ϵ = 56.0
# %%
outfile = "param_optim_polish_chi2full.txt"

p = [ϵ]
x₀ = [θ₀, Δθ₀, βᵣ₀]
lb = [0.9θ₀, 0.9Δθ₀, 0.9βᵣ₀]
ub = [1.1θ₀, 1.1Δθ₀, 1.1βᵣ₀]
prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼)
sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
x₀ = sol.u
χ² = sol.minimum
θ₀ = x₀[1]
Δθ₀ = x₀[2]
βᵣ₀ = x₀[3]
χ²s = stream.chi2_stream(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ, ic, r☼)
χ²c = stream.chi2_core(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ)
println(f, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5, "  ", χ², "  ",  χ²s, "  ", χ²c)
println("χ²:F-S-C:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5,
        "  ", χ², "  ",  χ²s, "  ", χ²c)




