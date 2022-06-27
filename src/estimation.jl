# Does a single estimation, starting from true parameter values estimated from
# the data.
function estimate_model(; vint::String = "0", only_likelihoods::Bool = false,
 			θ_init::Vector{T} = [0.32135,	5.6459,	14.855,	1.1614,  0.6486,   1.9196,	 14.771,
 					     -2.4895,	1,	    0.44004, 0.32415, 0.87235,  0.5,	 0.25921,
				  	     -9.1217,  80.267,	-13.647, 1.7296,  8.8188,   0.92622, 4.283,
				  	     4.9097,	0,      7.8609,	 7.739,	  0.011111, 6.8386,  6.5027,
				  	     0.028621]) where {T<:Float64}

	vars = matread("$path/data/DataToRun_pop09_boot.mat")

	## Data Renaming and Setup Pre-Bootstrap
	γ0vec     = vcat(quantile.(Gamma(0.5, 20), 0.005:0.01:0.895), 28:2:60, 64:4:100)
	δ_vec     = vcat(exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)), 3:2:20)

	data12    = vars["data12nopop"]
	data09    = vars["data09nopop"]
	bp        = vars["bpnopop"]

	distpara0 = vec(vars["distpara0"])

	x0 = θ_init[7:20]
	objectivefun(x) = objective(x, θ_init[7:20], distpara0, γ0vec, δ_vec, data12, data09, bp)

	x00 = deepcopy(x0)

	# TODO: investigate why these columns are being deleted
	x0 = x0[[1:2; 4:6; 8:end]]# [3, 7])
    #x0[[3,7]] .= 0.0
	if only_likelihoods
    	return get_likelihoods(x0, θ_init[7:20], distpara0, γ0vec, δ_vec, data12, data09, bp)
	end

	out =  objectivefun(x0)

	res     = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5), iterations = Int(1e5)))
	x, fval = res.minimizer, res.minimum

	θ = [x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
	CSV.write("bootstrap_estimates_$(vint).csv", Tables.table(θ))

	return θ
end

function get_likelihoods(x, x0, distpara0, γ0vec, δvec, data12, data09, bp)
	if x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0
        println("Bad news: Parameters out of bounds.")
        return Inf
    end
    xx = [x[1:2]; x0[3]; x[3:5]; x0[7]; x[6:(length(x0)-2)]]

	out = try
		full_model(xx, distpara0, γ0vec, δvec, data12, data09, bp)
	catch err
		@show typeof(err)
		Inf, Inf, Inf, Inf, Inf, Inf
	end
	f, distpara, fother, fWF, f1, f2 = out
    return f, f1, f2
end

function objective(x, x0, distpara0, γ0vec, δvec, data12, data09, bp)
    if x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0
        println("Bad news: Parameters out of bounds.")
        return Inf
    end
    xx = [x[1:2]; x0[3]; x[3:5]; x0[7]; x[6:(length(x0)-2)]]
	out = try
		full_model(xx, distpara0, γ0vec, δvec, data12, data09, bp)[1]
	catch err
		@show typeof(err)
		Inf
	end
    return out
end
