# ***************************************************************************************
# Code for "Match Quality, Search, and the Internet Market for Used Books" (2021)
# by Glenn and Sara Fisher Ellison
# ...
# This Julia implementation (Aug 2021 - June 2022) is by Reca Sarfati (sarfati@mit.edu),
# which is based on MATLAB code from Masao Fukui, Hongkai Zhang, & [others]
#
# Users should CTRL-F for all TODO's in this script to adjust settings for personal use.
# ***************************************************************************************
using CSV, DataFrames, Dates, Distributions, FixedEffectModels, MAT, Optim
using Random, RegressionTables, Roots, SparseArrays, Statistics

# TODO: Adjust flags below for what you want to run.
test_functions = true
use_parallel   = false
run_bootstrap  = false

# TODO: If bootstrap, choose between available modes: 1 or 2 (see description)
run_mode = 1

# TODO: Specify script parameters
vint      = "2022-06-14"
N_workers = 72  # No. workers to request from cluster
N_bs      = 200 # No. bootstrap iterations

# TODO: Set estimation hyperparameters
# ϵ / rounderr
#

# Build output folders if don't exist
path = dirname(@__FILE__)
!isdir("$path/plots")  && run(`mkdir $path/plots/`)
!isdir("$path/tables") && run(`mkdir $path/tables/`)
!isdir("$path/output_data") && run(`mkdir $path/output_data/`)

# Load functions
include("helpers.jl")
include("full_model.jl")

# Acquire parallel workers (specific to MIT server; requires modification)
if use_parallel
    # TODO: Accomodate code for personal computing environment.
end

# Test output of functions
if test_functions
    include("../test/helpers.jl")
end

# Bootstrap
if run_bootstrap
    #include("bootstrap_distpara.jl")
end
