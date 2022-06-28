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

	INPUT = "$path/../data/input"
	## Data Renaming and Setup Pre-Bootstrap
	data12    = vars["data12nopop"]
	data09    = vars["data09nopop"]
	bp        = vars["bpnopop"]

	data = Dict{Symbol,Dict{Symbol,Union{Vector{<:Number},Int64}}}()
	for (d_k, dataset) in zip([:of_09, :on_09, :on_12], [bp, data09, data12])
		data[d_k] = Dict{Symbol,Union{Vector{<:Number},Int64}}()
		for k in ["cdindex","first","cdid","mktsize","popular"]
	           data[d_k][Symbol(k)] = vecI64(dataset[k])
	    end
		for k in ["numlist","localint","conditiondif","p","obsweight","condition"]
	           data[d_k][Symbol(k)] = vecF64(dataset[k])
	   	end
		data[d_k][:N] = Int(dataset["N"])
		data[d_k][:M] = Int(dataset["M"])
	end
	data[:on_09][:basecond]  = vecF64(data09["basecond"])
	data[:on_12][:basecond]  = vecF64(data12["basecond"])
	data[:on_09][:pdif]      = vecF64(data09["pdif"])
	data[:on_12][:pdif]      = vecF64(data12["pdif"])
	data[:on_12][:disappear] = vecF64(data12["disappear"])

	@load "$INPUT/data_to_run.jld2" data



	#distpara0 = vec(vars["distpara0"])
	@load "$INPUT/distpara0.jld2" distpara0
	# distpara0 is 6-element vector:
	# 0.3346253805385
  	# 6.9715573158185
 	# 21.1108893646429
  	# 1.67012950793581
  	# 1.4463101539879
  	# 2.10965105553215

	x0 = θ_init[7:20]
	objectivefun(x) = objective(x, θ_init[7:20], distpara0, data12, data09, bp; parallel=parallel)

	x00 = deepcopy(x0)

	# TODO: investigate why these columns are being deleted
	x0 = x0[[1:2; 4:6; 8:end]]

	if only_likelihoods
    	return get_likelihoods(x0, θ_init[7:20], distpara0, data12, data09, bp; parallel=parallel)
	end
	res = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5),
		     	   iterations = Int(1e5), show_trace = true, store_trace = true))
	x, fval = res.minimizer, res.minimum

	θ = [x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
	CSV.write("bootstrap_estimates_$(vint).csv", Tables.table(θ))

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
