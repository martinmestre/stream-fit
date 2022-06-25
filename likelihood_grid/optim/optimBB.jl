"""Optimization trials."""

# %%
using AlgebraOfGraphics, CairoMakie
using PyCall
using GalacticOptim
using GalacticBBO


# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")

# %%
function χ²Core(x, p)
    β = x[1]*1.e-5
    θ = p[1]
    Δθ = p[2]
    ϵ = p[3]
    return stream.chi2_core(θ, Δθ, β, ϵ)
end

function χ²Stream(x, p)
    θ = x[1]
    Δθ = x[2]
    β = p[1]*1.e-5
    ϵ = p[2]
    return stream.chi2_stream(θ, Δθ, β, ϵ, ic, r☼)
end

function χ²Full(x, p)
    θ = x[1]
    Δθ = x[2]
    β = x[3]*1.e-5
    ϵ = p[1]
    return stream.chi2_full(θ, Δθ, β, ϵ, ic, r☼)
end

# %%
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
  -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122
# x₀ = [3.621802786938911822e+01, 2.745899649967862999e+01, 1.25]
θ₀, Δθ₀ = 38.741377485620355, 29.0303636347919
βᵣ₀ =  0.0002877703539252742*1.e5

# %%
# E = vcat([56],  range(60, 300, 25))
E = range(184., 300., step=2.)
outfile = "solutions_stream_core_BB.txt"
for ϵ in E
    global  θ₀, Δθ₀, βᵣ₀, x₀
    open(outfile, "a") do f
        # Core  χ²
        # p = [θ₀, Δθ₀, ϵ]
        # prob = OptimizationProblem(χ²Core, βᵣ₀, p, lb=0.8*βᵣ₀, ub=1.2*βᵣ₀)
        # sol = solve(prob, NLopt.LN_PRAXIS(), abstol=0.5e-1, reltol=5.0e-2, maxiters=50)
        # βᵣ₀ = sol.u[1]
        # χ² = sol.minimum
        # x₀ = [θ₀, Δθ₀, βᵣ₀]
        # println("χ²Core:",  ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3]*1.e-5, "  ", χ²)

        # # Stream χ²
        # p = [βᵣ₀, ϵ]
        # x₀ = [θ₀, Δθ₀]
        # prob = OptimizationProblem(χ²Stream, x₀, p, ic=ic, r☼=r☼, lb=0.8*x₀, ub=1.2*x₀)
        # sol = solve(prob, NLopt.LN_PRAXIS(), abstol=0.5e-1, reltol=5.0e-2, maxiters=300)
        # x₀ = sol.u
        # χ² = sol.minimum
        # θ₀ = x₀[1]
        # Δθ₀ = x₀[2]
        # x₀ = [θ₀, Δθ₀, βᵣ₀]
        # println(f, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3]*1.e-5, "  ", χ²)
        # println("χ²Stream:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3]*1.e-5, "  ", χ²)

        # Full χ²
        p = [ϵ]
        x₀ = [θ₀, Δθ₀, βᵣ₀]
        lb = [0.9θ₀, 0.9Δθ₀, 0.8βᵣ₀]
        ub = [1.1θ₀, 1.1Δθ₀, 1.5βᵣ₀]
        prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
        sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=1000)
        x₀ = sol.u
        χ² = sol.minimum
        θ₀ = x₀[1]
        Δθ₀ = x₀[2]
        βᵣ₀ = x₀[3]
        x₀ = [θ₀, Δθ₀, βᵣ₀]
        println(f, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3]*1.e-5, "  ", χ²)
        println("χ²Full:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3]*1.e-5, "  ", χ²)
    end
end # the file f is automatically closed after this block finishes


# # One step with NLopt full χ²
# prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼)
# sol = solve(prob, NLopt.LN_PRAXIS())

# # %%
# # Second step with NLopt full χ²
# x₀ = sol.u
# p = 65
# prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼)
# sol = solve(prob, NLopt.LN_PRAXIS())

# # %%
# # Third step with NLopt full χ²
# x₀ = sol.u
# p = 80.0
# prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼)
# sol = solve(prob, NLopt.LN_PRAXIS())

# # One step with core χ²
# x₀ = [1.25]
# p = [3.621802786938911822e+01, 2.745899649967862999e+01, 60.00]
# # prob = OptimizationProblem(χ²Core, x₀, p, ic=ic, r☼=r☼)
# # sol = solve(prob, NelderMead())

# # %%
# # One step with stream χ²
# x₀ = [p[1], p[2]]
# p= [1.5002414703369147, p[3]]
# prob = OptimizationProblem(χ²Stream, x₀, p, ic=ic, r☼=r☼)
# sol = solve(prob, NelderMead())

# # %%
# # One step with full constraint
# x₀ = [sol.u[1], sol.u[2], p[1]]
# p = [p[2]]
# print(x₀,  p)
# # %%
# prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼)
# sol = solve(prob, NelderMead())