"""Plot density profiles"""

using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
# using DelimitedFiles
# using CSV
# using DataFrames, DataFramesMeta


"""Parameters and initial conditions"""

param_file = "param_optim_sequentially_fermionmass_polished.txt"
mat = readdlm(param_file)
ϵ=mat[:,1]
θ=mat[:,2]
ω=mat[:,3]
β=mat[:,4]
W = θ+ω

ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122

f
pot_list = stream.pot_model(ϵ, θ, W, β)