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
#     display_name: Julia 1.7.2
#     language: julia
#     name: julia-1.7
# ---

# %%
"""Plot solutions from the optimization loop."""


using AlgebraOfGraphics, CairoMakie
using PyCall
using CSV
using DataFrames
using Showoff

# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")
potentials = pyimport("potential_classes")


# %%
# Load file of selected solutions.
df = DataFrame(CSV.File("solutions_stream_core_polished_with_GR_NLoptNelderMead.txt"; delim=" ", ignorerepeated=true))
println(df)


ϵ = [56, 100, 200, 300, 360]
println(ϵ)

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
      χ²[i] = df.chi2full[ind][1]
      param[i] = [ϵ[i], θ[i], W[i], β[i]]
      i += 1
end

# %%
# Compute MEPP selected solutions.
rar = [potentials.RAR(param[i]) for i = 1:length(ϵ)]
println("MEPP solutions computed.")


# %%
# Initial conditions needed for future use.
# Taking the initial conditions from the Galpy fit with fixed MW2014 potential.
ic_galpy = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
# Taking initial condition from an optimization with fixed potential ("fit_orbit_..._fixedpot.py")
ic = [148.87671997, 36.34168516, 7.95627538, -6.87147041, -12.48727587, -16.05002458]
r☼ = 8.122


# %%
# Core radii and masses; compute only if needed.
r_core = Vector{Float64}(undef, length(ϵ))
m_core = Vector{Float64}(undef, length(ϵ))
χ²stream = Vector{Float64}(undef, length(ϵ))

println("i ::: β ::: r_core ::: m_core ::: χ²stream")
for i = 1:length(ϵ)
      halo = potentials.RAR(param[i])
      temp = stream.get_core_GR(halo)
      r_core[i] = temp[1]
      m_core[i] = temp[2][1]
      χ²stream[i] = stream.chi2_stream(θ[i], Δθ[i], β[i], ϵ[i], ic, r☼)
      println(i,"  ",  ϵ[i],  "  ", r_core[i],"  ", m_core[i], "  ", χ²stream[i])
end

# %%
# Black hole constants needed.
G = 4.300923924e-6
c = 2.997925e5  # km s^-1
M_bh = 4.075e6  # M_sun
rₛ = 2.0*G*M_bh / c^2

# %%
# Make ρ(r) plot for the selected solutions.
update_theme!(fontsize=37)
plt = Vector{Layer}(undef, length(ϵ))
labels = (e->"$e").(ϵ)
style = [:dash, :dot, :dashdot, :dashdotdot, nothing]
println(labels)
length(plt)

# %%
for i = 1:length(ϵ)
      r = rar[i].r_s
      ρ_f = rar[i].mass_spline.derivative(1)
      ρ=ρ_f(r)
      df_plot = (x=r, y=ρ./(4π*r.^2),
                 Energy=fill(labels[i],length(r)),
                 )
      plt[i]=
            data(df_plot)*
            mapping(:x,:y;linestyle=:Energy=> L"$m_{\mathrm{fermion}}$ [keV/c²]")*
            visual(Lines, linewidth=2)
end

full_plt = sum(plt)

fig = Figure()
ag = draw!(fig[1,1], full_plt,
           axis=(xlabel=L"$r$ [kpc]",
                 ylabel=L"$ρ$ [$M_⊙$/kpc³]",
                 limits=((1.e-10,nothing),(1., nothing)),
                 xscale=log10, yscale=log10))

leg = AlgebraOfGraphics.compute_legend(ag)
println(fig[1,1])
ax = fig[1,1] |> contents |> first |> contents |> first # get the axis
axislegend(ax, leg...)

display(fig)

# %%
# Other way of doing the same plot:
size_inches = (4*2, 3*2)
size_pt = 72 .* size_inches
fig = Figure(resolution = size_pt, fontsize = 25)
gridpos = fig[1, 1]

f = draw!(gridpos, full_plt, axis=(xlabel=L"$r$ [kpc]",
ylabel=L"$ρ$ [$M_⊙$/kpc³]",
limits=((1.e-10,nothing),(1., nothing)),
xscale=log10, yscale=log10))
legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))

display(fig)
save("density_profiles.pdf", fig, pt_per_unit = 1)
println("ρ(r) plot done.")


# %%
# Suggestion from Ian Weaver (Julia slack)
r = rar[2].r_s

ϵs = ϵ
is = 1:5
labels = ["$x" for x ∈ ϵs]

function ρ_f(r, i)
      temp = rar[i].mass_spline.derivative(1)
      return temp(r)[1]/(4π*r^2)
end
# Each column is a separate profile
ρ_fs = ρ_f.(r, is')

# %%
# Build the DataFrame
df = DataFrame([r;; ρ_fs], [:r; [Symbol("ϵ_$ϵ") for ϵ ∈ ϵs]])

let

      size_inches = (4.2*2, 3*2)
      size_pt = 72 .* size_inches
      fig = Figure(resolution = size_pt, fontsize = 24)
      gridpos = fig[1, 1]
      grp = dims(1) => renamer(labels) => L"$m_{\mathrm{fermion}}$ [keV/c²]"
      plt = data(df) *
          mapping(:r => L"$r$ [kpc]", 2:ncol(df) .=> L"$ρ$ [$M_⊙$/kpc³]";
              color = grp,
              linestyle = grp
          ) *
          visual(Lines, linewidth=3)
      exp_rng=range(-9,1,step=2)
      xtickpos = (e->10.0.^e).(exp_rng)
      xticknames=replace.(Showoff.showoff(10. .^(exp_rng), :scientific),"1.0×"=> "" )
      f = draw!(gridpos, plt, axis=(
            limits=((1.e-10, 10^2),(1, nothing)),
            xscale=log10, yscale=log10, xgridvisible=false, ygridvisible=false,
            xticks = (xtickpos, xticknames)))
      legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))

      display(fig)
      save("density_profiles.pdf", fig, pt_per_unit = 1)
      println("ρ(r) plot done.")
end



# %%
# Plot in observable space.
i = 1
pot_list = stream.pot_model(ϵ[i], θ[i], W[i], β[i])
temp = stream.orbit_model(ic[1], ic[2], ic[3], ic[4], ic[5],ic[6], pot_list, r☼)
(ϕ₁, ϕ₂, d☼, μ_ra, μ_dec, v☼, x, y, z, v_circ) = temp

# %%
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

# %%
lab_sol = fill("MEPP model",length(ϕ₁))
lab_obs = fill("Observed", length(ϕ₁ₒ))
df_sol = DataFrame([ϕ₁, ϕ₂, d☼, μ_ra, μ_dec, v☼, lab_sol], [:ϕ₁, :ϕ₂, :d☼, :μ_ra, :μ_dec, :v☼, :lab_sol])
df_obs = DataFrame([ϕ₁ₒ, ϕ₂ₒ, ϕ₂ₛ, ϕ₂ᵢ, d☼ₒ, d☼ₛ, d☼ᵢ,
                    μ_raₒ, μ_raₛ, μ_raᵢ, μ_decₒ, μ_decₛ, μ_decᵢ,
                    v☼ₒ, v☼ₛ, v☼ᵢ, lab_obs],
                    [:ϕ₁ₒ, :ϕ₂ₒ, :ϕ₂ₛ, :ϕ₂ᵢ, :d☼ₒ, :d☼ₛ, :d☼ᵢ,
                    :μ_raₒ, :μ_raₛ, :μ_raᵢ, :μ_decₒ, :μ_decₛ, :μ_decᵢ,
                    :v☼ₒ, :v☼ₛ, :v☼ᵢ, :lab_obs])

# %%
let

      fig = Figure(resolution = size_pt, fontsize = 25)
      gridpos = fig[1, 1]
      plt_sol = data(df_sol)*
            mapping(:ϕ₁,:ϕ₂,color=:lab_sol=>"")*
            visual(Scatter, markersize=3)

      plt_obs = data(df_obs)*
            (mapping(:ϕ₁ₒ,:ϕ₂ₒ,color=:lab_obs=>"")+ mapping(:ϕ₁ₒ,:ϕ₂ₛ,color=:lab_obs=>"")+ mapping(:ϕ₁ₒ,:ϕ₂ᵢ,color=:lab_obs=>""))*
            visual(Lines, color="blue", linewidth=2, linestyle=:dot)

      f = draw!(gridpos, plt_sol+plt_obs, axis=(xlabel=L"$ϕ_1$ [degrees]",
                  ylabel=L"$ϕ_2$ [degrees]",
                  limits=((-90,10),(-4, 1))))
      legend!(gridpos, f; tellwidth=false, halign=:center, valign=:bottom, margin=(10, 10, 10, 10))


      display(fig)
      save("observables_position.pdf", fig, pt_per_unit = 1)
      println("plot done.")
end
