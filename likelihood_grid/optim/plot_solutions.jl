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

# %%
using AlgebraOfGraphics, CairoMakie
using PyCall
using CSV
using DataFrames

# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")
potentials = pyimport("potential_classes")


# %%
df = DataFrame(CSV.File("solutions_stream_core_merged.csv"; delim="  "))
# df = DataFrame(CSV.File("solutions_stream_core_GR_NLoptNelderMead.txt"; delim="  "))

# %%
df.ener_f = parse.(Float64, df.ener_f)
df.theta_0 = parse.(Float64, df.theta_0)
df.d_theta = parse.(Float64, df.d_theta)
df.beta_0 = parse.(Float64, df.beta_0)
println(df)

# %%
set_aog_theme!()
update_theme!(fontsize=45)
pdf = (x=df.ener_f, y=df.theta_0, z=df.d_theta, w= df.beta_0*1.e4, v=df.chi2)
plt = data(pdf)*visual(Lines, markersize=2)*(
      mapping(:x,:y)+mapping(:x,:z)*visual(color="red")+
      mapping(:x,:w)*visual(color="blue")+
      mapping(:x,:v)*visual(color="green"))
draw(plt, axis=(xlabel=L"ϵ [keV]", ylabel="Parameter",
                limits=(nothing,nothing))) |> display

# %%
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
      χ²[i] = df.chi2[ind][1]
      param[i] = [ϵ[i], θ[i], W[i], β[i]]
      i += 1
end
println(length(param))
# %%
println(param)
rar = [potentials.RAR(param[i]) for i = 1:length(ϵ)]

# %%
r_core = Vector{Float64}(undef, length(ϵ))
m_core = Vector{Float64}(undef, length(ϵ))
χ²stream = Vector{Float64}(undef, length(ϵ))

# Taking the initial conditions from the Galpy fit with fixed MW2014 potential.
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
# Taking initial condition from an optimization with fixed potential ("fit_orbit_..._fixedpot.py")
ic = [148.87671997, 36.34168516, 7.95627538, -6.87147041, -12.48727587, -16.05002458]
r☼ = 8.122

println("i ::: β ::: r_core ::: m_core ::: χ²stream")
for i = 1:length(ϵ)
      halo = potentials.RAR(param[i])
      temp = stream.get_core_GR(halo.r, halo.dnu_wrap, halo.mass_spline)
      r_core[i] = temp[1]
      m_core[i] = temp[2][1]
      χ²stream[i] = stream.chi2_stream(θ[i], Δθ[i], β[i], ϵ[i], ic, r☼)
      println(i,"  ",  ϵ[i],  "  ", r_core[i],"  ", m_core[i], "  ", χ²stream[i])
end

# %%
G = 4.300923924e-6
c = 2.997925e5  # km s^-1
M_bh = 4.075e6  # M_sun
rₛ = 2.0*G*M_bh / c^2

# %%
update_theme!(fontsize=37)
plt = Vector{Layer}(undef, length(ϵ))
g = ["black", "red", "blue", "green", "cyan"]
labels = (e->"ϵ=$e").(ϵ)
println(labels)
length(plt)

# %%
for i = 1:length(ϵ)
      ρ_f = rar[i].mass_spline.derivative(1)
      r = rar[i].r
      ρ=ρ_f(r)
      df_plot = (x=r/rₛ, y=ρ./(4π*r.^2), Energy=fill(labels[i],length(r)))
      plt[i]= data(df_plot)*mapping(:x,:y,color=:Energy=> L"$m_{\mathrm{fermion}}$ [keV/c²]")*visual(Lines, linewidth=4)
end

full_plt = sum(plt)

fig = Figure()
ag = draw!(fig[1,1], full_plt, axis=(xlabel=L"r [r_{\mathrm{Sch}}]",
            ylabel=L"$ρ$ [$M_⊙$/kpc³]",
            limits=((1.e-1,nothing),(1., nothing)),
            xscale=log10, yscale=log10))

leg = AlgebraOfGraphics.compute_legend(ag)
println(fig[1,1])
ax = fig[1,1] |> contents |> first |> contents |> first # get the axis
axislegend(ax, leg...)

display(fig)
println("done")

# %%
# Plot core mass for fixed values of θ₀ and W₀
# ϵ, θ₀, Δθ₀, β₀ = 380.0, 43.86771278429453, 29.572455417272877, 0.003519778013001537
#ϵ, θ₀, Δθ₀, β₀ = 377.9, 43.86728376941964, 29.57245716023231, 0.0035195067675640795
ϵ, θ₀, Δθ₀, β₀ = 370.0, 43.86633110156308, 29.572455135477718, 0.0035189127649460296
Δβ = 0.5β₀
W₀ = θ₀+Δθ₀
B = range(β₀+0.0Δβ, β₀+1.0Δβ, 50)
r_core = Vector{Float64}(undef, length(B))
m_core = Vector{Float64}(undef, length(B))
χ² = Vector{Float64}(undef, length(B))
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
    -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122
println("i   :::   β   :::   r_core   :::   m_core   :::   χ²")

for i = 1:length(B)
      β = B[i]
      # param = [ϵ, θ₀, W₀, β]
      pot_list = stream.pot_model(ϵ, θ₀, W₀, β)
      halo = pot_list[4]
      temp = stream.get_core_GR(halo.r, halo.dnu_wrap, halo.mass_spline)
      r_core[i] = temp[1]
      m_core[i] = temp[2][1]
      χ²[i] = stream.chi2_stream_potlist(pot_list, ic, r☼)
      println(i,"  ", β,  "  ", r_core[i],"  ", m_core[i])
end

# %%
update_theme!(fontsize=32)
ϵ, θ₀, Δθ₀ = 370.0, 43.86633110156308, 29.572455135477718
W₀ = θ₀+Δθ₀
β₀ = [0.0035189127649460296, 0.00352, 0.003761596403907824]
labels = (x->"$x").(β₀)
color = ["red", "cyan", "green"]
plt = Vector{Layer}(undef, length(β₀))
for i in range(1,length(β₀))
      param = [ϵ, θ₀, W₀, β₀[i]]
      halo = potentials.RAR(param)
      temp = stream.get_core_GR(halo.r, halo.dnu_wrap, halo.mass_spline)
      r_core = temp[1]
      m_core = temp[2][1]
      println(β₀[i],  "  ", r_core,"  ", m_core)
      ρ_f = halo.mass_spline.derivative(1)
      r = halo.r
      ρ=ρ_f(r)
      df = (x=r/rₛ, y=ρ./(4π*r.^2), β=fill(labels[i],length(r)) )
      plt[i] = data(df)*mapping(:x,:y, color=:β=>"β₀ (fixed θ₀,W₀,ϵ=370 keV)")*visual(Lines, linewidth=4)
end
full_plt = sum(plt)


fig = Figure()
ag = draw!(fig[1,1], full_plt, axis=(xlabel=L"r [r_{\mathrm{Sch}}]",
            ylabel=L"$ρ$ [$M_⊙$/kpc³]",
            limits=((1.e-1,nothing),(1., nothing)),
            xscale=log10, yscale=log10))

leg = AlgebraOfGraphics.compute_legend(ag)
ax = fig[1,1] |> contents |> first |> contents |> first # get the axis
axislegend(ax, leg...)

display(fig)

# %%
fig = Figure()
gridpos = fig[1, 1]

f = draw!(gridpos, full_plt, axis=(xlabel=L"r [r_{\mathrm{Sch}}]",
ylabel=L"$ρ$ [$M_⊙$/kpc³]",
limits=((1.e-1,nothing),(1., nothing)),
xscale=log10, yscale=log10))
legend!(gridpos, f; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10))

display(fig)

# %%
"""O-V mass."""
mₑ = 9.1093837015e-28/1.98847e33  # M⊙
Eₑ = 510.99895  # keV
E2m = mₑ/Eₑ
m2E = 1.0/E2m
function massOV(mₚ, m_f)
      return 0.384*mₚ^3/m_f^2
end
Eₚ= 1.22093e25
mₚ= Eₚ*E2m
E_f = 360.0
m_f = E_f*E2m
ov_mass = massOV(mₚ, m_f)
println(ov_mass)

function solve_fermion(ov_mass, mₚ)
      return sqrt(0.384mₚ^3/ov_mass)
end
ov_mass = 4.1e6
m_f = solve_fermion(ov_mass, mₚ)
E_f = m_f*m2E
println(E_f)