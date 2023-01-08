## TODO: Specify script parameters
vint    = "2023-01-08"
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

## TODO: Option to specify starting parameters
#θ_init = OrderedDict()
θ_init = OrderedDict(
    #=1=#	:α          => 14.7544,	#1  α 				[15.77 (1.28)]
    #=2=#	:Δ_p_out    => -2.5128,	#2  Δ_p_out 		[0.16 (0.02)]
    #=3=#	:γ_ns_shape => 1.0,		#3  γ_ns_shape *** FIXED
    #=4=#	:γ_ns_on_09 => 0.4247,	#4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
    #=5=#	:γ_ns_on_12 => 0.3275,	#5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
    #=6=#	:η          => 0.8302,	#6  η - 1 			[1.87 (0.08)]
    #=7=#	:r          => 0.5,     #7  r * 10     *** FIXED
    #=8=#	:R_p        => 0.2738,	#8  R_p 			[0.26 (0.02)]
    #=9=#	:c          => -9.0273,	#9  c * 10 			[-0.91 (0.07)]*10
    #=10=#	:γ_s_pop    => 78.7921,	#10 γ_s_pop * 100 	[0.80 (0.09)]
    #=11=#	:γ_ns_pop   => -14.592,	#11 γ_ns_pop * 10 	[-1.36 (0.19)]
    #=12=#	:s_R        => 1.7247,	#12 s_R 			[1.73 (0.09)]
    #=13=#	:μ_R        => 8.8787,  #13 μ_R / s_R		[15.25 (0.80) / 1.73 (0.09)]
    #=14=#	:R_q        => 0.9253)

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
bs_inds    = 1:2   # No. bootstrap iterations
seed       = true  # For replicating output / catching bugs
read_draws = ""

global path = dirname(@__FILE__)
include("$path/launch_script.jl")
