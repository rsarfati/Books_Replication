## TODO: Specify script parameters
vint    = "2023-01-08"
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
run_tests    = false # Test code matches MATLAB (for developers)
write_output = false # Saves output to file
estimation   = true # Estimate model
WFcal	     = false # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = true # Does NOT optimize; evaluates likelihood for given parameters

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

## TODO: Option to specify starting parameters
#θ_init = OrderedDict()
θ_init = OrderedDict(
    #=1=#	:α          => 15.4337,	#1  α 				[15.77 (1.28)]
    #=2=#	:Δ_p_out    => -2.684,	#2  Δ_p_out 		[0.16 (0.02)]
    #=3=#	:γ_ns_shape => 1.0,		#3  γ_ns_shape *** FIXED
    #=4=#	:γ_ns_on_09 => 0.408859,#4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
    #=5=#	:γ_ns_on_12 => 0.13914,	#5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
    #=6=#	:η          => 0.880357,#6  η - 1 			[1.87 (0.08)]
    #=7=#	:r          => 0.5,     #7  r * 10     *** FIXED
    #=8=#	:R_p        => 0.264626,#8  R_p 			[0.26 (0.02)]
    #=9=#	:c          => -8.7612,	#9  c * 10 			[-0.91 (0.07)]*10
    #=10=#	:γ_s_pop    => 78.7921,	#10 γ_s_pop * 100 	[0.80 (0.09)]
    #=11=#	:γ_ns_pop   => -14.592,	#11 γ_ns_pop * 10 	[-1.36 (0.19)]
    #=12=#	:s_R        => 1.84503,	#12 s_R 			[1.73 (0.09)]
    #=13=#	:μ_R        => 8.7488,  #13 μ_R / s_R		[15.25 (0.80) / 1.73 (0.09)]
    #=14=#	:R_q        => 0.92372)

global path = dirname(@__FILE__)
include("$path/launch_script.jl")
