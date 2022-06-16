using CSV, DataFrames, Dates, Distributed, Distributions, FixedEffectModels, MAT, Optim
using Random, RegressionTables, Roots, SparseArrays, Statistics

# Build output folders if don't exist
path = dirname(@__FILE__)
!isdir("$path/plots")       && run(`mkdir $path/plots/`)
!isdir("$path/tables")      && run(`mkdir $path/tables/`)
!isdir("$path/output_data") && run(`mkdir $path/output_data/`)

# ***************************************************************************************
# Code for "Match Quality, Search, and the Internet Market for Used Books" (2021)
# by Glenn and Sara Fisher Ellison
# ...
# This Julia implementation (Aug 2021 - June 2022) is by Reca Sarfati (sarfati@mit.edu),
# which is based on MATLAB code from Masao Fukui, Hongkai Zhang, & [others]
#
# Users should CTRL-F for all TODO's in this script to adjust settings for personal use.
# ***************************************************************************************

## TODO: Specify script parameters
vint    = "2022-06-15"
n_procs = 20  # No. workers to request from cluster
n_bs    = 200 # No. bootstrap iterations

## TODO: Adjust flags below for what you want to run.

test_functions = false  # Test code matches MATLAB (for developers)
parallel       = true  # Distribute work across multiple processes?
output_lik     = true  # Do you want to simply fetch the likelihood of a set of parameters?
estimation     = true  # Estimate model
run_bootstrap  = false # Run bootstrap for SEs?
run_mode       = 1     # Running bootstrap? Choose between modes 1 or 2 (see explanation).

## TODO: Set hyperparameters
rounderr = 0.025

# Add worker processes, load necessary packages on said workers
if parallel
    #addprocs(n_procs)
    @everywhere using CSV, DataFrames, Dates, Distributions, FixedEffectModels, MAT
    @everywhere using Optim, Random, RegressionTables, Roots, SparseArrays, Statistics
    @everywhere path = dirname(@__FILE__)
    println("Added worker processes!")
end

# Load functions
@everywhere include("$path/helpers.jl")
@everywhere include("$path/full_model.jl")

# Test function output, if you've been modifying code
if test_functions
    include("$path/../test/helpers.jl")
end

if output_lik
    include("$path/estimation.jl")
    f, f1, f2 = estimate_model(only_likelihoods = true)
    @save "likelihoods.jld2" f f1 f2
    #CSV.write("likelihoods.csv", )
end

# Estimate model from known parameters
if estimation
    include("$path/estimation.jl")
    estimate_model()
end

# Bootstrap
if run_bootstrap
    #include("bootstrap_distpara.jl")
end

#if parallel
#    rmprocs(workers())
#end
