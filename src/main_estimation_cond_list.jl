## TODO: Specify script parameters
vint    = "2023-01-22"
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
write_output = true # Saves output to file
estimation   = true # Estimate model
WFcal	     = true # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = true # Does NOT optimize; evaluates likelihood for given parameters

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

# Test code matches MATLAB (for developers)
#run_tests = false

## TODO: Option to specify starting parameters
#θ_init = OrderedDict()
#θ = Vector(CSV.read("../output/data/estimation_theta_cond_list_2023-01-10.jld2",
#                    DataFrame)[:,1])
# TODO: it is really weird that alpha is so precisely -1.
θ = [-0.999974033038977, -2.2943580777422588, 1.0, 0.12236310694992512,
     0.06694644933721287, 0.5771387813800685, 0.5, 0.829647280017562,
     -21.813036021385, 80.17094160788938, -17.579435604415647, 2.168397240948148,
     7.034165435127134, 0.9264674894656546, -0.0020925524094474544, 1.126346953052436,
     0.21681289559569902]
θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :α_c, :η_c,
                      :min_p] .=> θ)

# θ_init = OrderedDict(:α          	 =>	14.9942	,
# :Δ_p_out    	 =>	-2.4965	,
# :γ_ns_shape 	 =>	1.0	,
# :γ_ns_on_09 	 =>	0.08205	,
# :γ_ns_on_12 	 =>	0.03003	,
# :η          	 =>	0.57303	,
# :r          	 =>	0.5	,
# :R_p        	 =>	0.87539	,
# :c          	 =>	-9.00915	,
# :γ_s_pop    	 =>	78.7921	,
# :γ_ns_pop   	 =>	-14.592	,
# :s_R        	 =>	2.1005	,
# :μ_R        	 =>	6.83318	,
# :R_q        	 =>	0.9374	,
# :α_c  	 =>	0.08364	,
# :η_c         	 =>	3.9682	,
# :has_min_p	 =>	0.18424)


global path = dirname(@__FILE__)
include("$path/launch_script.jl")
