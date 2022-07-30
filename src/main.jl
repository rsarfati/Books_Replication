using CSV, DataFrames, Dates, Distributed, Distributions, FileIO
using FixedEffectModels, JLD2, MAT, Optim, OrderedCollections
using Printf, Random, Roots, SparseArrays, Statistics, UnPack
println("Packages loaded!")

# Build output folders if don't exist
global path   = dirname(@__FILE__)
global OUTPUT = "$path/../output/data"
global INPUT  = "$path/../input"

!isdir("$path/../output/")	&& run(`mkdir $path/../output/`)
!isdir("$OUTPUT/")			&& run(`mkdir $OUTPUT/`)
!isdir("$OUTPUT/../plots")	&& run(`mkdir $OUTPUT/../plots/`)
!isdir("$OUTPUT/../tables")	&& run(`mkdir $OUTPUT/../tables/`)

## TODO: Specify script parameters
vint    = "2022-07-29"
N_procs = 100	# No. workers to request from cluster
N_bs    = 50 	# No. bootstrap iterations

## TODO: Adjust flags below for what you want to run.
parallel   = true	# Distribute work across multiple processes?
run_tests  = false	# Test code matches MATLAB (for developers)
eval_only  = false	# Do you want to simply fetch the likelihood of a set of parameters?
estimation = true	# Estimate model
bootstrap  = false	# Run bootstrap for SEs?
run_mode   = :OPTIM	# Running bootstrap? Choose :OPTIM or :EVAL

# Add worker processes, load necessary packages on said workers
if parallel
	println("(1/2) Adding processes...")
    addprocs(N_procs)
    @everywhere using CSV, DataFrames, Dates, Distributed, Distributions, FileIO
    @everywhere using FixedEffectModels, JLD2, MAT, Optim, OrderedCollections
    @everywhere using Printf, Random, Roots, SparseArrays, Statistics, UnPack

    @everywhere global path = dirname(@__FILE__)
    @everywhere global OUTPUT = "$path/../output/data"
    @everywhere global INPUT  = "$path/../input"

    println("(2/2) Added $(length(workers())) worker processes!")
end

# Loads functions
@everywhere rounderr = 0.025
@everywhere include("$path/helpers.jl")
@everywhere include("$path/full_model.jl")
@everywhere include("$path/estimation.jl")
@everywhere include("$path/bootstrap.jl")

# Test function output (good idea if you've been modifying code)
if run_tests; include("$path/../test/helpers.jl") end

# Estimate model from known parameters
if estimation
	println("(1/2) Running estimation; eval_only = $(eval_only)")
    out = estimate_model(eval_only = eval_only, parallel = parallel)

    if eval_only # Only solve for likelihoods
        f, distpara, fother, fWF, f1, f2 = out
        @save "$OUTPUT/likelihoods_$(vint).jld2" f distpara fother fWF f1 f2
        @show f, f1, f2
    end
	println("(2/2) Finished running estimation.")
end

# Run bootstrap script
if bootstrap; run_bootstrap(N_bs = N_bs, parallel = parallel, eval_only = eval_only) end

# Release workers
if parallel; rmprocs(workers()); println("Workers released!") end
