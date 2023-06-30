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
"""Plot observables"""

using AlgebraOfGraphics, CairoMakie
using PyCall
using CSV
using DataFrames
using Showoff
using CubicSplines
using DelimitedFiles
include("bugfixMakie.jl")
# %%

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
# %%



# %%
# Parameters and initial conditions.
param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
θ, Δθ, β = readdlm(param_file)
ϵ = 56.0
W = θ+Δθ
param = [ϵ, θ, W, β]


ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122

orbit_nfw_file = "observable_orbit_NFW-MW.txt"
# %%


# Plot in observable space.

# GD-1 stream observables (Ibata polynomials).
ϕ₁ₒ = stream.Iba_phi1_np
ϕ₂ₒ = stream.Iba_phi2_np
d☼ₒ = stream.Iba_dhel_np
μ_raₒ = stream.Iba_mura_np
μ_decₒ = stream.Iba_mudec_np
v☼ₒ = stream.Iba_vhel_np

Δϕ₂ = fill(stream.sigma_array[1], length(ϕ₁ₒ))
Δd☼ = fill(stream.sigma_array[2], length(ϕ₁ₒ))
Δμ_ra = fill(stream.sigma_array[3], length(ϕ₁ₒ))
Δμ_dec = fill(stream.sigma_array[4], length(ϕ₁ₒ))
Δv☼ = fill(stream.sigma_array[5], length(ϕ₁ₒ))


ϕ₂ₛ = ϕ₂ₒ + Δϕ₂
ϕ₂ᵢ = ϕ₂ₒ - Δϕ₂
d☼ₛ = d☼ₒ + Δd☼
d☼ᵢ = d☼ₒ - Δd☼
μ_raₛ = μ_raₒ + Δμ_ra
μ_raᵢ = μ_raₒ - Δμ_ra
μ_decₛ = μ_decₒ + Δμ_dec
μ_decᵢ = μ_decₒ - Δμ_dec
v☼ₛ = v☼ₒ + Δv☼
v☼ᵢ = v☼ₒ - Δv☼

# MEPP solution.
pot_list = stream.pot_model(ϵ, θ, W, β)
temp = stream.orbit_model(ic..., pot_list, r☼)
# Cubic spline for solution:
boolϕ₁ = (temp[1].>-95 .&& temp[1].<15)
fϕ₂, fd☼, fμ_ra, fμ_dec, fv☼ = [CubicSpline(temp[1][boolϕ₁], temp[i][boolϕ₁]) for i=2:6]
ϕ₂ = fϕ₂(ϕ₁ₒ)
d☼ = fd☼(ϕ₁ₒ)
μ_ra = fμ_ra(ϕ₁ₒ)
μ_dec = fμ_dec(ϕ₁ₒ)
v☼ = fv☼(ϕ₁ₒ)

# Malhan (MWPotential2014) solution.
temp = readdlm(orbit_nfw_file)
boolϕ₁ = (temp[1,:].>-95 .&& temp[1,:].<15)
fϕ₂, fd☼, fμ_ra, fμ_dec, fv☼ = [CubicSpline(temp[1,:][boolϕ₁], temp[i,:][boolϕ₁]) for i=2:6]
ϕ₂ₘ = fϕ₂(ϕ₁ₒ)
d☼ₘ = fd☼(ϕ₁ₒ)
μ_raₘ = fμ_ra(ϕ₁ₒ)
μ_decₘ = fμ_dec(ϕ₁ₒ)
v☼ₘ = fv☼(ϕ₁ₒ)

df_obsmod = DataFrame([ϕ₁ₒ, ϕ₂ₒ, ϕ₂ₛ, ϕ₂ᵢ, d☼ₒ, d☼ₛ, d☼ᵢ,
                    μ_raₒ, μ_raₛ, μ_raᵢ, μ_decₒ, μ_decₛ, μ_decᵢ,
                    v☼ₒ, v☼ₛ, v☼ᵢ,
                    ϕ₂, d☼, μ_ra, μ_dec, v☼,
                    ϕ₂ₘ, d☼ₘ, μ_raₘ, μ_decₘ, v☼ₘ],
                    [:ϕ₁ₒ, :ϕ₂ₒ, :ϕ₂ₛ, :ϕ₂ᵢ, :d☼ₒ, :d☼ₛ, :d☼ᵢ,
                    :μ_raₒ, :μ_raₛ, :μ_raᵢ, :μ_decₒ, :μ_decₛ, :μ_decᵢ,
                    :v☼ₒ, :v☼ₛ, :v☼ᵢ,
                    :ϕ₂, :d☼, :μ_ra, :μ_dec, :v☼,
                    :ϕ₂ₘ, :d☼ₘ, :μ_raₘ, :μ_decₘ, :v☼ₘ])
# %%
labels = ["Observed±σ", "NFW-MW", "Fermionic-MW" ]
lw = 4
# %%

# Sky position plot
let
      size_inches = (6.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => ""
      plt = data(df_obsmod) *
          mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [2, 22, 17] .=> L"ϕ_2~[°]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=lw)
      plt2 = data(df_obsmod) *
      mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [3,4] .=> L"ϕ_2~[°]";
      ) *
      visual(Lines, linewidth=lw)
      println("plt=",plt+plt2)
      f = draw!(gridpos, plt+plt2, axis=(;limits=((-90,10),(-4, 1)),
            xgridvisible=false, ygridvisible=false))
      legend!(gridpos, f; tellwidth=false, halign=:center, valign=:bottom, margin=(10, 10, 10, 10), patchsize=(50,35))

      # Lines re-styling
      lineas = fig.content[1].scene.plots
      lineas[1].color = "black"
      lineas[4].linewidth = 2
      lineas[5].linewidth = 2
      lineas[2].linestyle = :dash
      lineas[3].linestyle = :dot

      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:4]
      for l in _lines
            l.linewidth = 4
      end
      _lines[1].color = "black"
      _lines[1].linestyle = :solid
      _lines[2].linestyle = :dash
      _lines[3].linestyle = :dot

      display(fig)
      save("paper_plots/observables_position.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end

# %%

# Proper motion (RA) plot
let
      size_inches = (6.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => ""
      plt = data(df_obsmod) *
          mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [8, 24, 19] .=> L"\tilde{μ}_{α}~[\mathrm{mas}~\mathrm{yr}^{-1}]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=lw)
      plt2 = data(df_obsmod) *
      mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [9,10] .=> L"\tilde{μ}_{α}~[\mathrm{mas}~\mathrm{yr}^{-1}]";
      ) *
      visual(Lines, linewidth=lw)
      f = draw!(gridpos, plt+plt2, axis=(;limits=((-90,10),(nothing, nothing)),
            xgridvisible=false, ygridvisible=false))
      legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(50,35))

      # Lines re-styling
      lineas = fig.content[1].scene.plots
      lineas[1].color = "black"
      lineas[4].linewidth = 2
      lineas[5].linewidth = 2
      lineas[2].linestyle = :dash
      lineas[3].linestyle = :dot

      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:4]
      for l in _lines
            l.linewidth = 4
      end
      _lines[1].color = "black"
      _lines[1].linestyle = :solid
      _lines[2].linestyle = :dash
      _lines[3].linestyle = :dot
      display(fig)
      save("paper_plots/observables_pmra.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end

# %%

# Proper motion (DEC) plot
let
      size_inches = (6.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => ""
      plt = data(df_obsmod) *
          mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [11, 25, 20] .=> L"μ_δ~[\mathrm{mas}~\mathrm{yr}^{-1}]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=lw)
      plt2 = data(df_obsmod) *
      mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [12,13] .=> L"μ_δ~[\mathrm{mas}~\mathrm{yr}^{-1}]";
      ) *
      visual(Lines, linewidth=lw)
      f = draw!(gridpos, plt+plt2, axis=(;limits=((-90,10),(nothing, nothing)),
            xgridvisible=false, ygridvisible=false))
      legend!(gridpos, f; tellwidth=false, halign=:center, valign=:top, margin=(10, 10, 10, 10), patchsize=(50,35))

      # Lines re-styling
      lineas = fig.content[1].scene.plots
      lineas[1].color = "black"
      lineas[4].linewidth = 2
      lineas[5].linewidth = 2
      lineas[2].linestyle = :dash
      lineas[3].linestyle = :dot

      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:4]
      for l in _lines
            l.linewidth = 4
      end
      _lines[1].color = "black"
      _lines[1].linestyle = :solid
      _lines[2].linestyle = :dash
      _lines[3].linestyle = :dot

      display(fig)
      save("paper_plots/observables_pmdec.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end

# %%
# Heliocentric distance plot
let
      size_inches = (6.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => ""
      plt = data(df_obsmod) *
          mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [5, 23, 18] .=> L"D~[\mathrm{kpc}]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=lw)
      plt2 = data(df_obsmod) *
      mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [6,7] .=> L"D~[\mathrm{kpc}]";
      ) *
      visual(Lines, linewidth=lw)
      f = draw!(gridpos, plt+plt2, axis=(;limits=((-90,10),(nothing, nothing)),
            xgridvisible=false, ygridvisible=false))
      legend!(gridpos, f; tellwidth=false, halign=:center, valign=:top, margin=(10, 10, 10, 10), patchsize=(50,35))

      # Lines re-styling
      lineas = fig.content[1].scene.plots
      lineas[1].color = "black"
      lineas[4].linewidth = 2
      lineas[5].linewidth = 2
      lineas[2].linestyle = :dash
      lineas[3].linestyle = :dot

      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:4]
      for l in _lines
            l.linewidth = 4
      end
      _lines[1].color = "black"
      _lines[1].linestyle = :solid
      _lines[2].linestyle = :dash
      _lines[3].linestyle = :dot

      display(fig)
      save("paper_plots/observables_heldist.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end
# %%

# Heliocentric radial velocity plot
let
      size_inches = (6.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => ""
      plt = data(df_obsmod) *
          mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [15, 26, 21] .=> L"v_h~[\mathrm{km~s^{-1}}]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=lw)
      plt2 = data(df_obsmod) *
      mapping(:ϕ₁ₒ => L"ϕ_1~[°]", [16] .=> L"v_h~[\mathrm{km~s^{-1}}]";
      ) *
      visual(Lines, linewidth=lw)
      f = draw!(gridpos, plt+plt2, axis=(;limits=((-90,10),(nothing, nothing)),
            xgridvisible=false, ygridvisible=false))

      legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(50,35))

       # Lines re-styling
      lineas = fig.content[1].scene.plots
      lineas[1].color = "black"
      lineas[1].linewidth = 2
      lineas[4].linewidth = 2
      lineas[2].linestyle = :dash
      lineas[3].linestyle = :dot

      leg = fig.content[2]
      _lines = leg.blockscene.children[1].plots[2:4]
      for l in _lines
            l.linewidth = 4
      end
      _lines[1].color = "black"
      _lines[1].linestyle = :solid
      _lines[2].linestyle = :dash
      _lines[3].linestyle = :dot

      display(fig)
      save("paper_plots/observables_helvel.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end