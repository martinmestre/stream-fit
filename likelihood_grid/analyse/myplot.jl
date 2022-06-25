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
"""Test PyCall.jl"""

# %%
using Pkg

Pkg.activate(".")
using PyCall
py"""
import sys
sys.path.insert(0, ".")
"""

# %%
plot_jl = pyimport("myplot")["plot_fun"]


# %%
dir = "serafin"
name = "likelihood_beta0_"
beta_array = collect(range(1.0e-5, 1.5e-5, 6))
beta_array[2] = 1.11e-05
println("beta_array = ", beta_array)
plot_jl(dir, name, beta_array)

