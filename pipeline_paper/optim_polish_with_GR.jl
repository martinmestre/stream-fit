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
df = DataFrame(CSV.File("solutions_stream_core_merged.csv"; delim="  "))

df.ener_f = parse.(Float64, df.ener_f)
df.theta_0 = parse.(Float64, df.theta_0)
df.d_theta = parse.(Float64, df.d_theta)
df.beta_0 = parse.(Float64, df.beta_0)
ϵ = [056, 100, 200, 300, 360]
println("ϵ = ", ϵ)

θ = Vector{Float64}(undef,length(ϵ))
Δθ = Vector{Float64}(undef,length(ϵ))
β = Vector{Float64}(undef,length(ϵ))
W = Vector{Float64}(undef,length(ϵ))
χ² = Vector{Float64}(undef,length(ϵ))
param = [Vector{Float64}(undef,4) for _ = 1:length(ϵ)]
i = 1
for E in ϵ
      global i
      ind = (df.ener_f.==E)
      θ[i] = df.theta_0[ind][1]
      Δθ[i] = df.d_theta[ind][1]
      W[i] = θ[i] + Δθ[i]
      β[i] = df.beta_0[ind][1]
      χ²[i] = df.chi2[ind][1]
      param[i] = [ϵ[i], θ[i], W[i], β[i]]
      i += 1
end
println(length(param))
println(param)
# %%
# Initial conditions from Galpy's MWPotential2014 model:
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122
# Taking initial condition from an optimization with fixed potential ("fit_orbit_..._fixedpot.py")
ic = [148.87671997, 36.34168516, 7.95627538, -6.87147041, -12.48727587, -16.05002458]
r☼ = 8.122

# %%
outfile = "solutions_stream_core_polished_with_GR_NLoptNelderMead.txt"
for i = 1:length(ϵ)
    open(outfile, "a") do f
        # Full χ²
        p = [ϵ[i]]
        x₀ = [θ[i], Δθ[i], β[i]*10^5]
        lb = [0.995θ[i], 0.995Δθ[i], 0.95β[i]*10^5]
        ub = [1.005θ[i], 1.005Δθ[i], 1.05β[i]*10^5]
        prob = OptimizationProblem(χ²Full, x₀, p, ic=ic, r☼=r☼, lb=lb, ub=ub)
        sol = solve(prob, NLopt.LN_NELDERMEAD(), reltol=5.0e-5)
        xₛ = sol.u
        χ² = sol.minimum
        θₛ = xₛ[1]
        Δθₛ = xₛ[2]
        βᵣₛ = xₛ[3]
        χ²s = stream.chi2_stream(θₛ, Δθₛ, βᵣₛ*1.e-5, ϵ[i], ic, r☼)
        χ²c = stream.chi2_core(θₛ, Δθₛ, βᵣₛ*1.e-5, ϵ[i])
        println(f, ϵ[i], "  ", xₛ[1], "  ", xₛ[2], "  ", xₛ[3] * 1.e-5, "  ", χ², "  ",  χ²s, "  ", χ²c)
        println("χ²:F-S-C:", ϵ[i], "  ", xₛ[1], "  ", xₛ[2], "  ", xₛ[3] * 1.e-5,
                "  ", χ², "  ",  χ²s, "  ", χ²c)

    end
end # the file f is automatically closed after this block finishes



