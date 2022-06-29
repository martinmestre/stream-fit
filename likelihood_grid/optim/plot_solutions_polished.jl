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

# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")
potentials = pyimport("potential_classes")


# %%
# Load file of selected solutions.
df = DataFrame(CSV.File("solutions_stream_core_polished_with_GR_NLoptNelderMead.txt"; delim=" ", ignorerepeated=true))
println(df)


ϵ = [056, 100, 200, 300, 360]
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
prinln("MEPP solutions comptued.")


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
labels = (e->"ϵ=$e").(ϵ)
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
# Plot in observable space.
