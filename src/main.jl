# ***************************************************************************************
# Code for "Match Quality, Search, and the Internet Market for Used Books" (2021)
# by Glenn and Sara Fisher Ellison
# ...
# This Julia implementation (Aug 2021 - June 2022) is by Reca Sarfati (sarfati@mit.edu),
# which is based on MATLAB code from Masao Fukui, Hongkai Zhang, & [others]
#
# Users should CTRL-F for all TODO's in this script to adjust settings for personal use.
# ***************************************************************************************
using CSV, DataFrames, Dates, Distributed, Distributions, FixedEffectModels, MAT, Optim
using Random, RegressionTables, Roots, SparseArrays, Statistics

## TODO: Specify script parameters

vint    = "2022-06-14"
n_procs = 2  # No. workers to request from cluster
n_bs    = 200 # No. bootstrap iterations

## TODO: Adjust flags below for what you want to run.

test_functions = true  # Test code matches MATLAB (for developers)
parallel       = true  # Distribute work across multiple processes?
estmation      = true  # Estimate model
run_bootstrap  = false # Run bootstrap for SEs?
run_mode = 1 # If running bootstrap, choose between modes 1 or 2 (see expltn.)

## TODO: Adjust estimation hyperparameters
# ϵ / rounderr

# Add worker processes, load necessary packages on said workers
if parallel
    addprocs(n_procs)
    @everywhere using CSV, DataFrames, Dates, Distributions, FixedEffectModels, MAT
    @everywhere using Optim, Random, RegressionTables, Roots, SparseArrays, Statistics
end

# Build output folders if don't exist
path = dirname(@__FILE__)
!isdir("$path/plots")  && run(`mkdir $path/plots/`)
!isdir("$path/tables") && run(`mkdir $path/tables/`)
!isdir("$path/output_data") && run(`mkdir $path/output_data/`)

# Load functions
@everywhere include("helpers.jl")
@everywhere include("full_model.jl")

# Test function output, if you've been modifying code
if test_functions
    include("../test/helpers.jl")
end

# Bootstrap
if run_bootstrap
    #include("bootstrap_distpara.jl")
end
