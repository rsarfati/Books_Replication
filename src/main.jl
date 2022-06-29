using CSV, DataFrames, Dates, Distributed, Distributions, FileIO, FixedEffectModels
using JLD2, MAT, Optim, OrderedCollections, Printf, Random, RegressionTables, Roots, SparseArrays, Statistics

# Build output folders if don't exist
global path   = dirname(@__FILE__)
global OUTPUT = "$path/../output"
global INPUT  = "$path/../input"

!isdir("$OUTPUT")        && run(`mkdir $OUTPUT`)
!isdir("$OUTPUT/plots")  && run(`mkdir $path/plots/`)
!isdir("$OUTPUT/tables") && run(`mkdir $path/tables/`)
!isdir("$OUTOUT/data")   && run(`mkdir $path/data/`)

## TODO: Specify script parameters
vint     = "2022-06-27"
n_procs  = 100  # No. workers to request from cluster
N_bs     = 200  # No. bootstrap iterations

## TODO: Adjust flags below for what you want to run.
parallel      = false  # Distribute work across multiple processes?
run_tests     = false  # Test code matches MATLAB (for developers)
eval_only     = true  # Do you want to simply fetch the likelihood of a set of parameters?
estimation    = false  # Estimate model
run_bootstrap = false  # Run bootstrap for SEs?
run_mode      = :OPTIM # Running bootstrap? Choose :OPTIM or :EVAL

# Add worker processes, load necessary packages on said workers
if parallel
    addprocs(n_procs)
    @everywhere using CSV, DataFrames, Dates, Distributed, Distributions
    @everywhere using FileIO, FixedEffectModels, JLD2, MAT, Optim, OrderedCollections, Printf
    @everywhere using Random, RegressionTables, Roots, SparseArrays, Statistics
    @everywhere global path = dirname(@__FILE__)
    @everywhere global OUTPUT = "$path/../output"
    @everywhere global INPUT  = "$path/../input"
    println("Added $(length(workers())) worker processes!")
end

# Loads functions
@everywhere rounderr = 0.025
@everywhere include("$path/helpers.jl")
@everywhere include("$path/full_model.jl")
@everywhere include("$path/estimation.jl")

# Test function output (good idea if you've been modifying code)
run_tests && include("$path/../test/helpers.jl")

# Only solve for likelihoods
if eval_only
    f, distpara, fother, fWF, f1, f2 = estimate_model(eval_only = true, parallel = parallel)
    @save "$OUTPUT/data/likelihoods_$(vint).jld2" f, distpara, fother, fWF, f1, f2
    @show f, f1, f2
end

# Estimate model from known parameters
estimation    && estimate_model()

# Run bootstrap script
run_bootstrap && include("$path/bootstrap_distpara.jl")

# Release workers
parallel      && rmprocs(workers())
