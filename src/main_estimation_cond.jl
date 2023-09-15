global path = dirname(@__FILE__)
include("$path/load_packages.jl")

## TODO: Specify script parameters
vint    = "2023-09-14"
spec    = :condition # Options are :standard, :condition, :cond_list
N_procs = 30	     # No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
write_output = true # Saves output to file
estimation   = true # Estimate model
WFcal	     = false # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = false # Does NOT optimize; evaluates likelihood for given parameters
grid_search  = false # meh
grid_search ? eval_only = true : nothing
VERBOSE      = true
run_tests    = false

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

# Test code matches MATLAB (for developers)
#run_tests = false

## TODO: Option to specify starting parameters
#θ = Vector(CSV.read("$path/../output/data/estimation_theta_standard_2023-01-10.csv",
#                    DataFrame)[:,1])

#@load "../output/data/estimation_results_2023_01_03.jld2"

#θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
#                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :α_c, :η_c] .=> [θ; 0.0; 0.0])

#@load "../input/condition_start.jld2"
@load "../output/data/estimation_results_condition_2023-07-30.jld2"
θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :α_c, :η_c] .=> θ)

include("$path/launch_script.jl")


## Grid search code below:
# points = OrderedDict{Symbol,Any}(
# 	:α          	 =>	4.0:6.0:16.0,
# 	:Δ_p_out    	 =>	-2.3,
# 	:γ_ns_shape 	 =>	1.0	,
# 	:γ_ns_on_09 	 =>	0.01:0.25:0.51,#0.023399967	,
# 	:γ_ns_on_12 	 =>	0.01:0.1:0.21,#0.005582159	,
# 	:η          	 =>	0.7,#0.6:0.15:0.9	,
# 	:r          	 =>	0.5,
# 	:R_p        	 =>	0.6,#0.3:0.3:0.9,#0.867298655	,
# 	:c          	 =>	-15.0:4.0:-7.0,#-9.203152599	,
# 	:γ_s_pop    	 =>	60.0:30.0:120.0,#73.30833091	,
# 	:γ_ns_pop   	 =>	-15.0,#-18.0:3.0:-12.0,#-14.59211909	,
# 	:s_R        	 =>	1.5:0.75:3.0,#2.439671533	,
# 	:μ_R        	 =>	6.0:2.5:8.5,#5.83948042	,
# 	:R_q        	 =>	0.927080228	,#0.8:0.02:0.96,#
# 	:α_c  	 		 =>	-1.0:2.0:5.0,#-0.008696885	,
# 	:η_c         	 =>	-1.0:1.5:3.5)
#
# global grid = Iterators.product(values(points)...)
# length(grid)
#
# global top_ten_params = Vector{Vector{Float64}}(undef, 10)
# global top_ten_llh    = repeat([Inf], 10)
# global ind_max = 1
# global ind_max_val = Inf
# global c = 0
#
# for g in grid
#
# 	global c += 1
# 	θ_g = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
#                       :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :α_c, :η_c] .=> [g...])
#
# 	out = estimate_model(θ_init = θ_g, eval_only = eval_only, spec = spec, parallel = parallel,
# 						 write_output = write_output, vint = vint, WFcal = WFcal)
# 	llh_g = out[1]
#
# 	# if c <= 10
# 	# 	top_ten_params[c] = out[2]
# 	# 	top_ten_llh[c]    = out[1]
# 	# 	ind_max = argmax(top_ten_llh)
# 	# 	ind_max_val = top_ten_llh[ind_max]
# 	# end
#
# 	if llh_g <= ind_max_val
# 		global top_ten_llh[ind_max] = llh_g
# 		global top_ten_params[ind_max] = out[2]
# 		global ind_max = argmax(top_ten_llh)
# 		global ind_max_val = top_ten_llh[ind_max]
# 	end
#
# 	if mod(c, 50) == 0
# 		@save "top_ten.jld2" top_ten_llh top_ten_params
# 		println("Parameters:")
# 		println(top_ten_params)
# 		println("LLHs:")
# 		println(top_ten_llh)
# 	end
# end
## Grid search code ends
