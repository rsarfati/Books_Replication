using CSV, DataFrames, Dates, DelimitedFiles, Distributed, Distributions, FileIO
using FixedEffectModels, JLD2, MAT, Optim, OrderedCollections
using Printf, Random, Roots, SparseArrays, Statistics, UnPack
println("Packages loaded!")

## Build output folders if don't exist
global OUTPUT = "$path/../output/data"
global INPUT  = "$path/../input"

!isdir("$path/../output/")	&& run(`mkdir $path/../output/`)
!isdir("$OUTPUT/")			&& run(`mkdir $OUTPUT/`)
!isdir("$OUTPUT/../plots")	&& run(`mkdir $OUTPUT/../plots/`)
!isdir("$OUTPUT/../tables")	&& run(`mkdir $OUTPUT/../tables/`)

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
    out = estimate_model(θ_init = θ_init, eval_only = eval_only, spec = spec, parallel = parallel,
						 write_output = write_output, vint = vint, WFcal = WFcal)
	println("(2/2) Finished running $(string(spec)) estimation.")
end

## Run bootstrap script
if bootstrap
	println("Starting Bootstrap for $(string(spec))! Indices: $(bs_inds)")
	run_bootstrap(θ_init = θ_init, bs_inds = bs_inds, seed = seed, spec = spec, vint = vint,
				  parallel = parallel, eval_only = eval_only, read_draws = read_draws)
	println("Bootstrap for $(string(spec)) complete! Output saved with vint = ``$vint``")
end

## Make output; only runs if you've already run bootstrap script!
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
#if parallel; rmprocs(workers()); println("Workers released!") end
