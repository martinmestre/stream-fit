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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
"""Analyse Likelihood data."""

# %%
using AlgebraOfGraphics, CairoMakie
using HDF5
using GLM
using PyCall

# %%
pushfirst!(PyVector(pyimport("sys")."path"), "")
stream = pyimport("stream")

# %%
# Testing and Learning!
pwd()
set_aog_theme!()
like_data=h5open("many_beta0s/serafin/model_beta0_1.24e-05.hdf5", "r")
table = read(like_data, "models_10002")
close(like_data)

x = Float64.(table["r"])
y = Float64.(table["mass"])
df = (x=x, y=y)
plt = data(df)* mapping(:x,:y)*visual(Lines)
draw(plt, axis=(xscale=log,)) # plot as heatmap (the default)


# %%
# Defining Circular velocity function (Wrapper to Python function).

function rot_vel_mw( ϵ, θ, W, β, ℜ)
    pot_list = stream.pot_model(ϵ, θ, W, β)
    dim = length(ℜ)
    v  = Array{Float64}(undef, dim)
    ix = 0
    for r in ℜ
        ix += 1
        v[ix] = stream.rot_vel_mw(pot_list, r)
    end
    return v
end
ℜ = range(1.0e-2, 50, 500)
ϵ, θ, Δθ, β = 56.0, 3.617308965727325187e+01, 2.740902073407499984e+01, 1.25e-5
W = θ+Δθ
v₁ = rot_vel_mw(ϵ, θ, W, β, ℜ)
ϵ, θ, Δθ, β = 56.0, 3.628885544245152772e+01, 2.751120739351376443e+01, 1.25e-5
W = θ+Δθ
v₂ = rot_vel_mw(ϵ, θ, W, β, ℜ)

set_aog_theme!()
update_theme!(fontsize=35)
df = (x=ℜ, y=v₁, z=v₂)
plt = data(df)*visual(Lines, markersize=2)*(mapping(:x,:y)+mapping(:x,:z)*visual(color="red"))
draw(plt, axis=(xlabel=L"$r$ [kpc]", ylabel=L"$V_{\mathrm{rot}}$ [km/s]",
                limits=((0, 50), (0, 300)))) |> display

# %%
# Reading data to make the plot.

# like_data = h5open("many_beta0s/serafin/model_beta0_1.24e-05.hdf5", "r")
like_data = h5open("beta0_1.25e-5/serafin/model_beta0_1.25e-05.hdf5", "r")

dim = length(keys(like_data))
θ₀ = Array{Float64}(undef, dim)
Δθ₀ = Array{Float64}(undef, dim)
W₀ = Array{Float64}(undef, dim)
β₀ = Array{Float64}(undef, dim)
loglike = Array{Float64}(undef, dim)
mass_vir = Array{Float64}(undef, dim)
radio_vir = Array{Float64}(undef, dim)
v_circ_rar  = Array{Float64}(undef, dim)
i = 0
for key in keys(like_data)
    global g
    global i
    i += 1
    g = like_data[key]
    θ₀[i] = read_attribute(g,"theta_0")
    Δθ₀[i] = read_attribute(g, "d_theta")
    β₀[i] = read_attribute(g, "beta_0")
    loglike[i] = read_attribute(g, "loglike")
    radio_vir[i] = last(g["r"])
    mass_vir[i] = last(g["mass"])
    W₀[i] = θ₀[i] + Δθ₀[i]
end
print("If i=",i, "=10^6, then everything is OK.")

close(like_data)


# %%
# Looking at the Likelihood in the parameter 2D slice.
set_aog_theme!()
update_theme!(fontsize=30)
bool = loglike .> -50
df = (x=θ₀[bool], y=Δθ₀[bool], z=loglike[bool], w=mass_vir[bool]/10^11)
plt = data(df)* mapping(:x,:y; color=:z)*visual(Scatter, markersize=10)
draw(plt, axis=(xlabel=L"θ_0", ylabel=L"W_0-θ_0"))

# %%
# Compute the dark matter mass in order to be at 1σ from Table 3 at Gibbons+14 for r =50 kpc.

mass_mw = 2.9e11
σ = 0.4e11
unit_mass_Allen = 2.32e7
mass_disks = 2.0*1700.0*unit_mass_Allen
mass_bulge = 460.0*unit_mass_Allen
mass_barions = mass_disks+mass_bulge
mass_halo = mass_mw - σ - mass_barions
println("Disks mass = ", mass_disks/1.e11, " x 10^11 M_⊙")
println("Barions mass = ", mass_barions/1.e11, " x 10^11 M_⊙")
println("Halo mas = ", mass_halo/1.e11, " x 10^11 M_⊙")

# %%
# Looking at the Likelihood along a 1D slice.
bool = loglike .> -50
df = (x=θ₀[bool], y=Δθ₀[bool], z=loglike[bool], w=mass_vir[bool]/10^11)
plt = data(df)* mapping(:x,:w)*visual(Scatter, markersize=2)
draw(plt, axis=(xlabel=L"θ_0", ylabel=L"M_{v}")) |> display
minarg = argmin(abs.(mass_vir[bool]-mass_halo*ones(Float64,length(θ₀[bool]))))

# %%
# Linear Fit to the crest.
x = θ₀[bool]
y = Δθ₀[bool]
linearRegressor = lm(@formula(y~x), df)


# %%
# Plot the linear fit.
linearFit = predict(linearRegressor)
df_lf = (u=θ₀[bool], v=linearFit)
plt = data(df)* mapping(:x,:y; color=:z)*visual(Scatter)
plt_lf = data(df_lf)*mapping(:u,:v)*visual(Lines; color="red")
u = range(35, 37, length=100)
f(x)=0.781486+0.736136*x
h(x)=0.784292+0.736041*x # Use h for the analysis. # without 1/3 factor in Pression
h(x)=0.768237+0.736953*x
v₁ = f.(u)
v₂ = h.(u)
df_lf_yo = (u=u, v=v₁)
df_lf_me = (u=u, v=v₂)
plt_lf_yo = data(df_lf_yo)*mapping(:u,:v)*visual(Lines; color="cyan")
plt_lf_me = data(df_lf_me)*mapping(:u,:v)*visual(Lines; color="brown")
draw(plt+plt_lf_yo+plt_lf+plt_lf_me)


# %%
# Computation of χ² along the crest.
ic = [1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
      -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01]
r☼ = 8.122
h(x)=0.784292+0.736041*x
h(x)=0.768237+0.736953*x-0.001  # With 1/3 factor in Pression corrected.
θ = range(35, 37, 50)
Δθ = h.(θ)
W = θ + Δθ
β = 1.25e-5
ϵ = 56.0
dim = length(θ)
χs² = Array{Float64}(undef, dim)
for i in range(1, dim)
    χ²s[i], r, mass, nu = stream.chi2(θ[i], Δθ[i], β, ϵ, ic, r☼)
end


# %%
# Plot χ² along the crest.
set_aog_theme!()
update_theme!(fontsize=35)
df = (x=θ, y=χ², group=fill("crest",length(θ)))
dfs = (x=θ, y=χ²s, group=fill("shifted",length(θ)))
plt = (data(df)+data(dfs))*visual(Lines, markersize=2)*mapping(:x,:y,color=:group)
draw(plt, axis=(xlabel=L"θ", ylabel=L"χ²", limits=((35, 37),nothing))) |> display

# %%
