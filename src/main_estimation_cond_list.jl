## TODO: Specify script parameters
vint    = "2023-01-08"
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

## TODO: Option to specify starting parameters
#θ_init = OrderedDict()
θ_init = OrderedDict()

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

global path = dirname(@__FILE__)
include("$path/launch_script.jl")
