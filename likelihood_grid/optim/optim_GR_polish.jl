"""Optimization trials."""

# %%
using AlgebraOfGraphics, CairoMakie
using PyCall
using GalacticOptim
using GalacticNLopt

# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")

# %%
function χ²Stream(x, p)
    θ = x[1]
    Δθ = x[2]
    β = p[1] * 1.e-5
    ϵ = p[2]
    return stream.chi2_stream(θ, Δθ, β, ϵ, ic, r☼)
end

function χ²Full(x, p)
    θ = p[1]
    Δθ = p[2]
    β = x[1] * 1.e-5
    ϵ = p[3]
    return stream.chi2_full(θ, Δθ, β, ϵ, ic, r☼)
end



# %%
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122

θ₀, Δθ₀, βᵣ₀ = 43.86633110156308, 29.572455135477718, 0.0037615964039078246*1e5
E = [370.0]

# %%
outfile = "solutions_stream_core_GR_polish_NLoptNelderMead.txt"
for ϵ in E
    println("ϵ=",ϵ)
    global θ₀, Δθ₀, βᵣ₀
    open(outfile, "a") do f
        # Full χ²
        p = [θ₀, Δθ₀, ϵ]
        x₀ = [βᵣ₀]
        lb = [0.9βᵣ₀]
        ub = [1.1βᵣ₀]
        prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
        sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
        x₀ = sol.u
        χ² = sol.minimum
        βᵣ₀ = x₀[1]
        x₀ = [βᵣ₀]
        println(f, "χ²Full ::: ", ϵ, " ::: ", x₀[1] * 1.e-5, "  ", χ²)
        println("χ²Full ::: ", ϵ, " :::  ",  x₀[1] * 1.e-5, "  ", χ²)
    end
end # the file f is automatically closed after this block finishes

# %%
outfile = "solutions_stream_core_GR_polish_NLoptNelderMead.txt"
for ϵ in E
    println("ϵ=",ϵ)
    global θ₀, Δθ₀, βᵣ₀
    open(outfile, "a") do f
        # Stream χ²
        p = [βᵣ₀, ϵ]
        x₀ = [θ₀, Δθ₀]
        lb = [0.9θ₀, 0.9Δθ₀]
        ub = [1.1θ₀, 1.1Δθ₀]
        prob = OptimizationProblem(χ²Stream, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
        sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
        x₀ = sol.u
        χ² = sol.minimum
        θ₀ = x₀[1]
        Δθ₀ = x₀[2]
        x₀ = [θ₀, Δθ₀]
        println(f, ϵ, "  ", x₀[1], "  ", x₀[2], "  ", p[1]*1.e-5, "  ", χ²)
        println("χ²Stream:", ϵ, "  ", x₀[1], "  ", x₀[2], "  ", p[1]*1.e-5, "  ", χ²)
    end
end # the   file f is automatically closed after this block finishes



