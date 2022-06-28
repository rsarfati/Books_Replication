# Does a single estimation, starting from true parameter values estimated from
# the data.
function estimate_model(; vint::String = "0", only_likelihoods::Bool = false,
			parallel::Bool = true,
 			θ_init::Vector{T} = [0.32135, 5.6459, 14.855, 1.1614, 0.6486,
								 1.9196,
								 # 7:20
								 14.771,  #1  α
								 -2.4895, #2  Δ_p_out
								 1,       #3  γ_ns_shape     *** FIXED
								 0.44004, #4  -> xx[4] * 10  = (γ_ns_on_09 / γ_ns_shape) * 9.5 ^ (-η) / r
								 0.32415, #5  -> xx[5] * 10  = (γ_ns_on_12 / γ_ns_shape) * 8.0 ^ (-η) / r
								 0.87235, #6  η - 1          *** I am infering there is meant to be a -1 transformation
								 0.5,     #7  r / 0.1        *** FIXED
								 0.25921, #8  R_p
								 -9.1217, #9  c / 0.1
				  	       		 80.267,  #10 γ_s_pop / 0.01
								 -13.647, #11 γ_ns_pop / 0.1
								 1.7296,  #12 s_R
								 8.8188,  #13 μ_R / s_R
								 0.92622, #14 1 - R_q
								 4.283,   #15 γ_s_shape
								 4.9097,  #16 γ_s_on_09
								 0,       #17 γ_s_on_12
								 7.8609,  #18 σ_δ
								 7.739,   #19 γ_ns_of_09_std
								0.011111, #20 βlocal ( = γ_ns_of_09_loc / γ_ns_of_09_std)
								 6.8386, 6.5027, 0.028621]) where {T<:Float64}

	INPUT  = "$path/../data/input"
	OUTPUT = "$path/../data/output"

	# Load data
	@load "$INPUT/data_to_run.jld2" data
	@load "$INPUT/distpara0.jld2"   distpara0
	# distpara0 is 6-element vector:
	# 0.3346253805385
  	# 6.9715573158185
 	# 21.1108893646429
  	# 1.67012950793581
  	# 1.4463101539879
  	# 2.10965105553215

	x0 = θ_init[7:20]
	objectivefun(x) = objective(x, θ_init[7:20], distpara0, data[:on_12],
								data[:on_09], data[:of_09]; parallel=parallel)

	x00 = deepcopy(x0)

	# TODO: investigate why these columns are being deleted
	x0 = x0[[1:2; 4:6; 8:end]]

	if only_likelihoods
    	return get_likelihoods(x0, θ_init[7:20], distpara0, data[:on_12],
							   data[:on_09], data[:of_09]; parallel=parallel)
	end

	# Optimize the objective function
	res = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5),
		     	   iterations = Int(1e5), show_trace = true, store_trace = true))
	x, fval = res.minimizer, res.minimum

	# Splice up the estimates
	θ = [x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]] # Fix parameters 3,7

	# optimize: x0 ->
	# estimates = [i; x[1:2]; x00[3]; x[3:5]; x00[7]; x[6:(length(x00)-2)]]

	# boot_out  = [i, θ]
	# xx        = θ[:, 2:15]
	# distpara  = boot_out[:, 16:21]

	# est = [boot_out[:, 16:21];  # 1:6 [γ_s_shape, γ_s_on_09, γ_s_on_12, σ_δ, γ_ns_of_09_std, βlocal (=γ_ns_of_09_loc/γ_ns_of_09_std)]
	# 	   xx[1];                 # 7 α
	# 	   xx[2] / (1 + xx[1]);   # 8 Δ_p_out
	# 	   xx[3];                 # 9 (FIXED: γ_ns_shape)
	# 	   xx[4] * 10 * xx[7] / 10 / 9.5 ^ (-xx[6] - 1); # 10
	# 	   xx[5] * 10 * xx[7] / 10 / 8.0 ^ (-xx[6] - 1); # 11
	# 	   xx[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1];      # 12:17 (η, r, R_p, c, γ_s_pop, γ_ns_pop)
	# 	   0;         # 18
	# 	   0;         # 19
	# 	   xx[12:13]; # 20:21
	# 	   xx[14]]    # 22
	CSV.write("$OUTPUT/bootstrap_estimates_$(vint).csv", Tables.table(θ))

	return θ
end

function get_likelihoods(x, x0, distpara0, data12, data09, bp; parallel::Bool=true)
	if x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0
        println("Bad news: Parameters out of bounds.")
        return Inf
    end
    xx = [x[1:2]; x0[3]; x[3:5]; x0[7]; x[6:(length(x0)-2)]]
	out = try
		full_model(xx, distpara0, data12, data09, bp; parallel=parallel)
	catch err
		#print(err)
		throw(err)
		Inf, Inf, Inf, Inf, Inf, Inf
	end
	f, distpara, fother, fWF, f1, f2 = out
    return f, f1, f2
end

function objective(x, x0, distpara0, data12, data09, bp; parallel::Bool=true)
    if x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0
        println("Bad news: Parameters out of bounds.")
        return Inf
    end
    xx = [x[1:2]; x0[3]; x[3:5]; x0[7]; x[6:(length(x0)-2)]]
	out = try
		full_model(xx, distpara0, data12, data09, bp; parallel=parallel)[1]
	catch err
		#print(err)
		throw(err)
		Inf
	end
    return out
end
