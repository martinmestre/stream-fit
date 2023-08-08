# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
"""Plot density profiles"""

# %%
using Pkg
Pkg.activate(".")
using AlgebraOfGraphics, CairoMakie
using DelimitedFiles
using PyCall
using Showoff
using DataFrames, DataFramesMeta
include("bugfixMakie.jl")

# %%
# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)

# %%
"""Parameters and initial conditions"""

# %%
param_file = "param_optim_sequentially_fermionmass_polished.txt"
mat = readdlm(param_file)
ϵ=mat[:,1]
θ=mat[:,2]
ω=mat[:,3]
β=mat[:,4]
W = θ+ω

# %%
ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122

# %%
#pot_list = stream.pot_model(ϵ, θ, W, β)

# Core radii and masses; compute only if needed.
r_core = Vector{Float64}(undef, length(ϵ))
m_core = Vector{Float64}(undef, length(ϵ))
χ²stream = Vector{Float64}(undef, length(ϵ))
param = [Vector{Float64}(undef,4) for _ = 1:length(ϵ)]
println("i ::: β ::: r_core ::: m_core ::: χ²stream")
for i = eachindex(ϵ)
      param[i] = [ϵ[i], θ[i], W[i], β[i]]
    #   halo = potentials.RAR(param[i])
    #   temp = stream.get_core_GR(halo)
    #   r_core[i] = temp[1]
    #   m_core[i] = temp[2][1]
    #   χ²stream[i] = stream.chi2_stream(θ[i], ω[i], β[i], ϵ[i], ic, r☼)
    #   println(i,"  ",  ϵ[i],  "  ", r_core[i],"  ", m_core[i], "  ", χ²stream[i])
end
# %%

# Compute MEPP selected solutions.
rar = [potentials.RAR(param[i]) for i = 1:length(ϵ)]
println("MEPP solutions computed.")
# %%

# Black hole constants needed.
G = 4.300923924e-6  # kpc (km/s)^2 M_sun^(-1)
c = 2.997925e5  # km s^(-1)
M_bh = 4.075e6  # M_sun
rₛ = 2.0*G*M_bh / c^2
# %%

# Suggestion from Ian Weaver (Julia slack)
r = rar[2].r_s

ϵs = ϵ
is = 1:5

labels = ["$x" for x ∈ ϵs]

function ρ_f(r, i)
      return rar[i].rho_spl(r)[1]
end
# Each column is a separate profile
ρ_fs = ρ_f.(r, is')
# %%

# Build the DataFrame
df = DataFrame([r;; ρ_fs], [:r; [Symbol("ϵ_$ϵ") for ϵ ∈ ϵs]])

let
      set_aog_theme!()
      update_theme!(Axis=(topspinevisible=true, rightspinevisible=true,
      topspinecolor=:darkgray, rightspinecolor=:darkgray,))
      size_inches = (5.2*2, 3.5*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => L"$m$ [keV/c²]"
      plt = data(df) *
          mapping(:r => L"$r$ [kpc]", 2:ncol(df) .=> L"$ρ$ [$M_⊙$/kpc³]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=5)
      exp_rng=range(-9,1,step=2)
      xtickpos = (e->10.0.^e).(exp_rng)
      xticknames=replace.(Showoff.showoff(10. .^(exp_rng), :scientific),"1.0×"=> "" )
      f = draw!(gridpos, plt, axis=(
            limits=((1.e-10, 10^2),(1, nothing)),
            xscale=log10, yscale=log10, xgridvisible=false, ygridvisible=false,
            xticks = (xtickpos, xticknames), yticksmirrored = true))
      legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(50,35))
      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:5]
      for l in _lines
            l.linewidth = 4
      end
      exp_rng=range(0,10,step=2)
      xtickpos = (e->10.0.^e).(exp_rng)
      xticknames=replace.(Showoff.showoff(10. .^(exp_rng), :scientific),"1.0×"=> "" )
      ax2 = Axis(fig[1, 1], xticklabelcolor = :black, xaxisposition = :top,
      limits=((1.e-10/rₛ, 10^2/rₛ),(1, nothing)), xlabel=L"$r$ [Schwarzschild radius]",
      xscale=log10, yscale=log10, xgridvisible=false, ygridvisible=false,
      xticks = (xtickpos, xticknames))
      hidespines!(ax2)
      hideydecorations!(ax2)
      display(fig)
      save("paper_plots/density_profiles.pdf", fig, pt_per_unit = 1)
      println("ρ(r) plot done.")
end