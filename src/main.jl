using CSV, DataFrames, Dates, Distributed, Distributions, FileIO, FixedEffectModels
using JLD2, MAT, Optim, Printf, Random, RegressionTables, Roots, SparseArrays, Statistics

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
vint     = "2022-06-26"
n_procs  = 100   # No. workers to request from cluster
N_bs     = 200   # No. bootstrap iterations
rounderr = 0.025 # Round error for stores

## TODO: Adjust flags below for what you want to run.
parallel      = true   # Distribute work across multiple processes?
run_tests     = false  # Test code matches MATLAB (for developers)
output_lik    = true   # Do you want to simply fetch the likelihood of a set of parameters?
estimation    = false  # Estimate model
run_bootstrap = true   # Run bootstrap for SEs?
run_mode      = :OPTIM # Running bootstrap? Choose :OPTIM or :EVAL

# Add worker processes, load necessary packages on said workers
if parallel
    addprocs(n_procs)
    @everywhere using CSV, DataFrames, Dates, Distributed, Distributions
    @everywhere using FileIO, FixedEffectModels, JLD2, MAT, Optim, Printf
    @everywhere using Random, RegressionTables, Roots, SparseArrays, Statistics
    @everywhere path     = dirname(@__FILE__)
    @everywhere rounderr = 0.025
    println("Added $(length(workers())) worker processes!")
end

# Loadsfunctions
@everywhere include("$path/helpers.jl")
@everywhere include("$path/full_model.jl")
@everywhere include("$path/estimation.jl")

# Test function output (good idea if you've been modifying code)
run_tests ? include("$path/../test/helpers.jl") : nothing

# Only solve for likelihoods
if output_lik
    f, f1, f2 = estimate_model(only_likelihoods = true)
    @save "likelihoods.jld2" f f1 f2
end

# Estimate model from known parameters
estimation    ? estimate_model()                       : nothing
# Run bootstrap script
run_bootstrap ? include("$path/bootstrap_distpara.jl") : nothing
# Release workers
parallel      ? rmprocs(workers())                     : nothing
