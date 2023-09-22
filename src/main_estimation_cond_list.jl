## TODO: Specify script parameters
vint    = "2023-01-22-b" # give today's date.
spec    = :cond_list # Options are :standard, :condition, :cond_list
N_procs = 30	 	 # No. workers to request from cluster

# TODO: Adjust flags below for what you want to run.
parallel     = true # Distribute work across multiple processors
write_output = true # Saves output to file
estimation   = true # Estimate model
WFcal	       = true # Grab welfare statistics
bootstrap    = false # Run bootstrap for SEs
eval_only    = true # Does NOT optimize; evaluates likelihood for given parameters

# TODO: Bootstrap flags
bs_inds     = 1:2 # No. bootstrap iterations
seed        = true # For replicating output / catching bugs
make_output = false # Prints pretty tables from bootstrap
read_draws  = ""

# Test code matches MATLAB (for developers)
#run_tests = false

 θ = [-0.9999453672410827,
      -2.2942922232742897,
      1.0,
      0.12235876222593084,
      0.0669359955701015,
      0.5771252840463006,
      0.5,
      0.8296363291652983,
      -21.817618569571284,
      80.17470393944173,
      -17.578892711699847,
      2.168330411073163,
      7.034013805867226,
      0.9264716291025366,
      -0.002083207311053787,
      1.1259145492381246,
      0.2167855707297015]

θ_init = OrderedDict([:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
                      :R_p, :c , :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :α_c, :η_c,
                      :min_p] .=> θ)

global path = dirname(@__FILE__)
include("$path/launch_script.jl")
