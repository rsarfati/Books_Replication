iter = length(ARGS) == 0 ? "" : "_$(ARGS[1])"
serv = length(ARGS) == 0 ? "" : "_$(ARGS[2])"

global path = dirname(@__FILE__)
include("$path/load_packages.jl")

## TODO: Specify script parameters
vint    = "2023-09-21" * serv * iter
spec    = :standard # Options are :standard, :condition, :cond_list
N_procs = 10	 	# No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
write_output = true # Saves output to file
estimation   = true # Estimate model
WFcal	     = false # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = false # Does NOT optimize; evaluates likelihood for given parameters
VERBOSE      = false # whether or not to print statements
run_tests	 = false
out_to_log   = true

# Settings for optimization
max_iter = 200

# if out_to_log
# 	run(`touch $LOGS/$(vint)_log.csv`)
# end

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = false # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

# Test code matches MATLAB (for developers)
# run_tests = false

## TODO: Option to specify starting parameters
@load "../output/data/estimation_results_standard_2023-09-22.jld2"
θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q] .=> θ)

include("$path/launch_script.jl")
