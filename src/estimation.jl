# Does a single estimation, starting from true parameter values estimated from
# the data.
function estimate_model(; vint::String = "0",
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
	objectivefun(x) = objective(x, x0, distpara0, γ0vec, δ_vec, data12, data09, bp)

	x00 = x0

	# TODO: investigate why these columns are being deleted
	#deleteat!(x0, [3, 7])

	res     = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5), iterations = Int(1e5)))
	x, fval = res.minimizer, res.minimum

	θ = [i, x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
	CSV.write("bootstrap_estimates_$(vint).csv", Tables.table(θ))

	return θ
end

function objective(x, x0, distpara0, γ0vec, δvec, data12, data09, bp)
    if x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0
        return Inf
    end

    xx = vcat(x[1:2], x0[3], x[3:5], x0[7], x[6:(length(x0)-2)])

    @assert x[3] == xx[4]
    @assert x[4] == xx[5]
    @assert x[12] == xx[14]

    # if xx[4]<0 || xx[5]<0 || xx[14]>1 || xx[14]< 0
    #     return Inf
    # end
    return full_model(xx, distpara0, γ0vec, δvec, data12, data09, bp)
end
