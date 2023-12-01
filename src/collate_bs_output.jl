include("load_packages.jl")
parallel = true
N_procs  = 15
VERBOSE = true
write_output = true
WFcal	     = true # Grab welfare statistics
eval_only    = true # Does NOT optimize; evaluates likelihood for given parameters

day = "2023-09-21"

global vint    = day
spec    = :standard # Options are :standard, :condition, :cond_list
N_procs = 15	 	# No. workers to request from cluster

if parallel
    println("(1/2) Adding processes...")
    addprocs(N_procs)
    @everywhere using CSV, DataFrames, Dates, Distributed, Distributions, FileIO
    @everywhere using JLD2, MAT, Optim, OrderedCollections
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

N_bs = 100
N_θ  = 14
N_dp = 6

θ_bs       = zeros(Float64, N_bs, N_θ)
distpara_bs = zeros(Float64, N_bs, N_dp)
llh_bs     = zeros(Float64, N_bs)

# @load "$INPUT/data_to_run.jld2" data
# @load "$INPUT/bootstrap_indices.jld2" bootindex

θ_names = [:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
					  :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q]
distpara_names = [:γ_s_shape, :γ_s_on_09, :γ_s_on_12, :σ_δ, :γ_ns_of_09_std, :γ_ns_of_09_loc]

for i=1:N_bs
	@load "$OUTPUT/estimation_results_standard_2023-09-21_run=$i.jld2" θ_i llh_i distpara_i
	θ_init = OrderedDict(θ_names .=> θ_i)
	global vint = day *"_run=$i"

	@load "$INPUT/data_to_run.jld2" data
	@load "$INPUT/bootstrap_indices.jld2" bootindex

	data_i		 = deepcopy(data)
	bootindex_i  = deepcopy(bootindex)
	θ_i_t 		 = deepcopy(θ_i)
	distpara_i_t = deepcopy(distpara_i)

	# Transformations on θ
	θ_i_t[2] = -θ_i[2]/(1.0 + θ_i[1])
	θ_i_t[4] = θ_i[4] * 10.0 * θ_i[7] / 10.0 / 9.5 ^ (-θ_i[6] - 1.0)
	θ_i_t[5] = θ_i[5] * 10.0 * θ_i[7] / 10.0 / 8.0 ^ (-θ_i[6] - 1.0)
	θ_i_t[6] = θ_i[6] + 1.0
	θ_i_t[7] = θ_i[7] * 0.1
	θ_i_t[9] = θ_i[9] * 0.1
	θ_i_t[10] = θ_i[10] * 0.01
	θ_i_t[11] = θ_i[11] * 0.1
	θ_i_t[14] = 1.0 - θ_i[14]

	# Transformations on distpara
	distpara_i_t[1:3] .= abs.(distpara_i[1:3])
	distpara_i_t[5]    = abs(distpara_i[5])
	distpara_i_t[6]    = distpara_i_t[5] * distpara_i[6]

	θ_bs[i,:]		 = θ_i_t
	distpara_bs[i,:] = distpara_i_t
	llh_bs[i] = llh_i

	# Welfare computation
	out = estimate_model(θ_init = θ_init,
						 data = index_data(data_i, bootindex_i[i,:]),
					     distpara0 = distpara_i,
						 eval_only = true, spec = :standard, parallel = true,
						 write_output = true, vint = day * "_run=$i", WFcal = true,
						 VERBOSE = true)

	println("Transformed and saved welfare computation for i=$i.")
end

CSV.write("$OUTPUT/bootstrap_theta_$(string(spec))_$(day)_all.csv", Tables.table(θ_bs; header=string.(θ_names)))
CSV.write("$OUTPUT/bootstrap_distpara_$(string(spec))_$(day)_all.csv", Tables.table(distpara_bs;header=string.(distpara_names)))
@save     "$OUTPUT/bootstrap_full_estimation_results_$(string(spec))_$(day)_all.jld2" θ_bs llh_bs distpara_bs

## Quick re-run
# spec    = :standard # Options are :standard, :condition, :cond_list
# N_procs = 15	 	# No. workers to request from cluster
#
# parallel     = true # Distribute work across multiple processors
# write_output = true # Saves output to file
# estimation   = true # Estimate model
# WFcal	     = false # Grab welfare statistics
# bootstrap    = false # Run bootstrap for SEs
# eval_only    = false # Does NOT optimize; evaluates likelihood for given parameters
# VERBOSE      = true # whether or not to print statements
# run_tests	 = false
# out_to_log   = true
#
# @load "../output/data/estimation_results_2023_01_03.jld2"
# θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
#                       :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q] .=> θ_bs[argmin(llh_bs),:])
#
# @everywhere global path   = "/bbkinghome/sarfati/ra_work/Books_Replication/src"
# @everywhere global OUTPUT = "$path/../output/data"
# @everywhere global INPUT  = "$path/../input"
#
# if parallel
#   println("(1/2) Adding processes...")
#   addprocs(N_procs)
#   @everywhere using CSV, DataFrames, Dates, Distributed, Distributions, FileIO
#   @everywhere using JLD2, MAT, Optim, OrderedCollections
#   @everywhere using Printf, Random, Roots, SparseArrays, Statistics, UnPack
#
#   @everywhere global path   = "/bbkinghome/sarfati/ra_work/Books_Replication/src"
#   @everywhere global OUTPUT = "$path/../output/data"
#   @everywhere global INPUT  = "$path/../input"
#
#   println("(2/2) Added $(length(workers())) worker processes!")
# end
#
# ## Load functions on all processors
# @everywhere rounderr = 0.025
# @everywhere include("helpers.jl")
# @everywhere include("full_model.jl")
# @everywhere include("estimation.jl")
# @everywhere include("bootstrap.jl")
#
# estimate_model(θ_init = θ_init, eval_only = eval_only, spec = spec, parallel = parallel,
# 					 write_output = write_output, vint = vint, WFcal = WFcal, VERBOSE = VERBOSE)
