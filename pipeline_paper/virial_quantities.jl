"""Compute virial quantities of our solution."""

using Pkg
Pkg.activate(".")
using PyCall
using CSV
using DelimitedFiles

pushfirst!(PyVector(pyimport("sys")."path"), "")
importLib = pyimport("importlib")
stream = pyimport("stream")
potentials = pyimport("potential_classes")
importLib.reload(stream)
importLib.reload(potentials)

# Parameters and initial conditions.
param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
θ, Δθ, β = readdlm(param_file)
ϵ = 56.0
W = θ+Δθ
param = [ϵ, θ, W, β]


ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = readdlm(ic_file)
r☼ = 8.122
# %%

pot_list = stream.pot_model(ϵ, θ, W, β)
# %%

rfk = pot_list[4]
@show rfk.r_max==rfk.r_s[end]
r_vir = rfk.r_max
m_vir = rfk.mass_wrap(r_vir)[1]
@show r_vir m_vir;

mass_mw = 2.9e11
σ = 0.9e11
unit_mass_Allen = 2.32e7
mass_disks = 2.0*1700.0*unit_mass_Allen
mass_bulge = 460.0*unit_mass_Allen
mass_barions = mass_disks+mass_bulge
mass_halo = mass_mw - σ - mass_barions
println("Disks mass = ", mass_disks/1.e11, " x 10^11 M_⊙")
println("Barions mass = ", mass_barions/1.e11, " x 10^11 M_⊙")
println("Halo mas = ", mass_halo/1.e11, " x 10^11 M_⊙")