using CSV, DataFrames, Dates, DelimitedFiles, Distributed, Distributions, FileIO
using FixedEffectModels, JLD2, MAT, Optim, OrderedCollections
using Printf, Random, Roots, SparseArrays, Statistics, UnPack
println("Packages loaded!")

## Build output folders if don't exist
global path   = dirname(@__FILE__)
global OUTPUT = "$path/../output/data"
global INPUT  = "$path/../input"

!isdir("$path/../output/")	&& run(`mkdir $path/../output/`)
!isdir("$OUTPUT/")			&& run(`mkdir $OUTPUT/`)
!isdir("$OUTPUT/../plots")	&& run(`mkdir $OUTPUT/../plots/`)
!isdir("$OUTPUT/../tables")	&& run(`mkdir $OUTPUT/../tables/`)

## TODO: Specify script parameters
vint    = "2023-01-05"
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
run_tests    = false # Test code matches MATLAB (for developers)
write_output = true # Saves output to file
estimation   = true # Estimate model
WFcal	     = false # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = false # Does NOT optimize; evaluates likelihood for given parameters
make_output  = false

# TODO: Bootstrap flags
bs_inds = 1:2   # No. bootstrap iterations
seed    = true  # For replicating output / catching bugs

## Add worker processes, load necessary packages on said workers
if parallel
    println("(1/2) Adding processes...")
    addprocs(N_procs)
    @everywhere using CSV, DataFrames, Dates, Distributed, Distributions, FileIO
    @everywhere using FixedEffectModels, JLD2, MAT, Optim, OrderedCollections
    @everywhere using Printf, Random, Roots, SparseArrays, Statistics, UnPack

    @everywhere global path   = dirname(@__FILE__)
    @everywhere global OUTPUT = "$path/../output/data"
    @everywhere global INPUT  = "$path/../input"

    println("(2/2) Added $(length(workers())) worker processes!")
end

## Load functions on all processors
@everywhere rounderr = 0.025
@everywhere include("$path/helpers.jl")
@everywhere include("$path/full_model.jl")
@everywhere include("$path/estimation.jl")
@everywhere include("$path/bootstrap.jl")

## Test function output (good idea if you've been modifying code)
if run_tests; include("$path/../test/helpers.jl") end

## Estimate model, starting from known parameters
if estimation
    println("(1/2) Estimating $(string(spec)) model; eval_only = $(eval_only)")
    out = estimate_model(eval_only = eval_only, spec = spec, parallel = parallel,
						 write_output = write_output, vint = vint, WFcal = WFcal)
	println("(2/2) Finished running $(string(spec)) estimation.")
end

## Run bootstrap script
if bootstrap
	println("Starting Bootstrap for $(string(spec))! Indices: $(bs_inds)")
	run_bootstrap(bs_inds = bs_inds, seed = seed, spec = spec, vint = vint,
				  parallel = parallel, eval_only = eval_only)
	println("Bootstrap for $(string(spec)) complete! Output saved with vint = ``$vint``")
end

if make_output
	θ_bs = vcat([Vector(CSV.read("$OUTPUT/bs_llh_theta_$(string(spec))_$(vint)_run=$i.csv",
								 header=true, DataFrame)[2:end,1])' for i=bs_inds]...)
	distpara_bs = vcat([Vector(CSV.read("$OUTPUT/bs_distpara_$(string(spec))_$(vint)_run=$i.csv",
								 header=true, DataFrame)[:,1])' for i=bs_inds]...)
	# Compute and format results from estimation
	b_boot = output_statistics(boot = θ_bs, distpara = distpara_bs,
	 						   #boot_out = "$OUTPUT/bootstrap_welfare_$(vint).csv",
	                           vint = vint, write_out = true)[1]
	make_table_results(b_boot;
		table_title = "$OUTPUT/../tables/estimates_$(string(spec))_$(vint)_eval_only=$(eval_only).tex")
end

## Release workers
if parallel; rmprocs(workers()); println("Workers released!") end
