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
# Initial conditions from Galpy's MWPotential2014 model:
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122
# Taking initial condition from an optimization with fixed potential ("fit_orbit_..._fixedpot.py")
ic = [148.87671997, 36.34168516, 7.95627538, -6.87147041, -12.48727587, -16.05002458]
r☼ = 8.122
# x₀ = [3.621802786938911822e+01, 2.745899649967862999e+01, 1.25]
# E = vcat([56],  range(60, 300, 25))
# θ₀, Δθ₀, βᵣ₀ = 39.285705296552436, 29.261703198023056, 0.0004816356819529165 * 1.e5
# E = range(218.3, 400.0, step=0.1)
θ₀, Δθ₀, βᵣ₀ = 36.224542166107774, 27.46389908753495, 1.2472408164492195
E = range(56.0, 380.0, step=0.1)
println("ϵ ∈ E = ", E)
# %%
outfile = "solutions_stream_core_GR_NLoptNelderMead.txt"
for ϵ in E
    global θ₀, Δθ₀, βᵣ₀, x₀
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
        lb = [0.999θ₀, 0.99Δθ₀, 0.99βᵣ₀]
        ub = [1.01θ₀, 1.01Δθ₀, 1.1βᵣ₀]
        prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
        sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
        x₀ = sol.u
        χ² = sol.minimum
        θ₀ = x₀[1]
        Δθ₀ = x₀[2]
        βᵣ₀ = x₀[3]
        x₀ = [θ₀, Δθ₀, βᵣ₀]
        χ²s = stream.chi2_stream(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ, ic, r☼)
        χ²c = stream.chi2_core(θ₀, Δθ₀, βᵣ₀*1.e-5, ϵ)
        println(f, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5, "  ", χ², "  ",  χ²s, "  ", χ²c)
        println("χ²:F-S-C:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", x₀[3] * 1.e-5,
                "  ", χ², "  ",  χ²s, "  ", χ²c)

    end
end # the file f is automatically closed after this block finishes



