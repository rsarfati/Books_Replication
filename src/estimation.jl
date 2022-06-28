# Does a single estimation, starting from true parameter values estimated from
# the data.
function estimate_model(; vint::String = "0", only_likelihoods::Bool = false,
			parallel::Bool = true,
 			θ_init::Vector{T} = [0.32135, 5.6459, 14.855, 1.1614, 0.6486,
								 1.9196, 14.771, -2.4895,
								 1, 0.44004, #==#
								 0.32415, 0.87235, 0.5, 0.25921,#==# -9.1217,
				  	       		 80.267, -13.647, 1.7296, 8.8188, 0.92622,
								 4.283, 4.9097, 0, 7.8609, 7.739,
								 0.011111, 6.8386, 6.5027, 0.028621]) where {T<:Float64}

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
	θ = [x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
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
