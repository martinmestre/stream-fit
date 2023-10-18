"""Compare core density of my solution vs the one that fits S-stars."""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using DelimitedFiles
using PyCall
using Showoff
using DataFrames, DataFramesMeta
include("bugfixMakie.jl")
#%%

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)
#%%

matlab_file = "/home/mmestre/casa/work/2020/precession/rar-potential/tables_matlab/final/56/Mc_3500000/rho.txt"

mat = readdlm(matlab_file)
r=mat[:,1]/10^3
ρ=mat[:,2]*10^9
# %%

"""Parameters and initial conditions"""
m = 56
param_file = "sol_optim_pot_m$(Int(m)).txt"
θ, ω, β = vec(readdlm(param_file))
W = θ+ω
param = [m,θ,W,β]
rar = potentials.RAR(param)
#%%
function ρ_f(r)
    return rar.rho_spl(r)[1]
end
bool_r = r.<10^-4.
r = r[bool_r]
ρ_gd1 = ρ_f.(r)
ρ = ρ[bool_r]

Δρ=ρ-ρ_gd1
@show Δρ