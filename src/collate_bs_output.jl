include("load_packages.jl")
parallel = true
N_procs  = 15

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
dispara_bs = zeros(Float64, N_bs, N_dp)
llh_bs     = zeros(Float64, N_bs)

@load "$INPUT/data_to_run.jld2" data

for i=1#:N_bs
	@load "$OUTPUT/estimation_results_standard_2023-09-21_run=$i.jld2" θ_i llh_i distpara_i

	θ_i_t 		 = deepcopy(θ_i)
	distpara_i_t = deepcopy(distpara_i)

	# Transformations on θ
	θ_i_t[2] = -θ_i[1]/(1.0 + θ_i[2])
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

	# Welfare computation

	out = estimate_model(θ_init = θ_i,
						 data = index_data(data, bootindex[i,:]),
					     distpara0 = distpara_i,
						 eval_only = true, spec = :standard, parallel = true,
						 write_output = true, vint = "2023-09-21", WFcal = true, VERBOSE = true)
	println("Transformed and saved welfare computation for i=$i.")
end
@show θ_bs[1,:]
@show distpara_bs[1,:]
