## TODO: Specify script parameters
vint    = "2023-01-08"
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

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

## TODO: Option to specify starting parameters
#θ_init = OrderedDict()
θ_init = OrderedDict(:α          	 =>	14.9942	,
:Δ_p_out    	 =>	-2.4965	,
:γ_ns_shape 	 =>	1.0	,
:γ_ns_on_09 	 =>	0.08205	,
:γ_ns_on_12 	 =>	0.03003	,
:η          	 =>	0.57303	,
:r          	 =>	0.5	,
:R_p        	 =>	0.87539	,
:c          	 =>	-9.00915	,
:γ_s_pop    	 =>	78.7921	,
:γ_ns_pop   	 =>	-14.592	,
:s_R        	 =>	2.1005	,
:μ_R        	 =>	6.83318	,
:R_q        	 =>	0.9374	,
:α_c  	 =>	0.08364	,
:η_c         	 =>	3.9682	,
:has_min_p	 =>	0.18424)


global path = dirname(@__FILE__)
include("$path/launch_script.jl")
