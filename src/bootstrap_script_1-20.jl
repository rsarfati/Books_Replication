iter = length(ARGS) == 0 ? "" : "_$(ARGS[1])"
serv = length(ARGS) == 0 ? "" : "_$(ARGS[2])"

global path = dirname(@__FILE__)
include("$path/load_packages.jl")

## TODO: Specify script parameters
vint    = "2023-09-21" * serv * iter
spec    = :standard # Options are :standard, :condition, :cond_list
N_procs = 15	 	# No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
write_output = false # Saves output to file within estimation
estimation   = false # Estimate model
WFcal	     = false # Grab welfare statistics
bootstrap    = true # Run bootstrap for SEs
eval_only    = false # Does NOT optimize; evaluates likelihood for given parameters
VERBOSE      = true # whether or not to print statements
run_tests	 = false
out_to_log   = true

# Settings
max_iter = 200

# if out_to_log
# 	run(`touch $LOGS/$(vint)_log.csv`)
# end

# TODO: Bootstrap flags
bs_inds     = 1:20 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

# Test code matches MATLAB (for developers)
# run_tests = false

## TODO: Option to specify starting parameters
#θ = Vector(CSV.read("../output/data/estimation_results_2023_01_03.jld2",#estimation_theta_standard_2023-01-10.csv",
#                    DataFrame)[:,1])
@load "../output/data/estimation_results_2023_01_03.jld2"
θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q] .=> θ)
#θ_init = OrderedDict()
# θ_init = OrderedDict(
#     #=1=#	:α          => 15.4337,	#1  α 				[15.77 (1.28)]
#     #=2=#	:Δ_p_out    => -2.684,	#2  Δ_p_out 		[0.16 (0.02)]
#     #=3=#	:γ_ns_shape => 1.0,		#3  γ_ns_shape *** FIXED
#     #=4=#	:γ_ns_on_09 => 0.408859,#4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
#     #=5=#	:γ_ns_on_12 => 0.13914,	#5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
#     #=6=#	:η          => 0.880357,#6  η - 1 			[1.87 (0.08)]
#     #=7=#	:r          => 0.5,     #7  r * 10     *** FIXED
#     #=8=#	:R_p        => 0.264626,#8  R_p 			[0.26 (0.02)]
#     #=9=#	:c          => -8.7612,	#9  c * 10 			[-0.91 (0.07)]*10
#     #=10=#	:γ_s_pop    => 78.7921,	#10 γ_s_pop * 100 	[0.80 (0.09)]
#     #=11=#	:γ_ns_pop   => -14.592,	#11 γ_ns_pop * 10 	[-1.36 (0.19)]
#     #=12=#	:s_R        => 1.84503,	#12 s_R 			[1.73 (0.09)]
#     #=13=#	:μ_R        => 8.7488,  #13 μ_R / s_R		[15.25 (0.80) / 1.73 (0.09)]
#     #=14=#	:R_q        => 0.92372)

include("$path/launch_script.jl")
