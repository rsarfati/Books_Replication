# Does a single estimation, starting from true parameter values estimated from
# the data.
function estimate_model(; vint::String = "0", only_likelihoods::Bool = false,
			parallel::Bool = true,
 			θ_init::Vector{T} = [#[0.32135, 5.6459, 14.855, 1.1614, 0.6486,
								# 1.9196,
								 # 7:20
								 14.771,  #1  α
								 -2.4895, #2  Δ_p_out
								 1,       #3  γ_ns_shape     *** FIXED
								 0.44004, #4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
								 0.32415, #5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
								 0.87235, #6  η - 1          *** I am infering there is meant to be a -1 transformation
								 0.5,     #7  r * 10         *** FIXED
								 0.25921, #8  R_p
								 -9.1217, #9  c * 10
				  	       		 80.267,  #10 γ_s_pop * 100
								 -13.647, #11 γ_ns_pop * 10
								 1.7296,  #12 s_R
								 8.8188,  #13 μ_R / s_R
								 0.92622, #14 1 - R_q
								 4.283,   #15 γ_s_shape
								 4.9097,  #16 γ_s_on_09
								 0,       #17 γ_s_on_12
								 7.8609,  #18 σ_δ
								 7.739,   #19 γ_ns_of_09_std
								 0.011111 #20 βlocal ( = γ_ns_of_09_loc / γ_ns_of_09_std)
								 #6.8386,  #21
								 #6.5027,  #22
								 #0.028621 #23
								 ],
			θ_fix = OrderedDict{Symbol,Float64}(:r => 0.5, :γ_ns_shape => 1)) where {T<:Float64}

	# File path constants
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

	# Ordered list of symbols
	θ_sym = [:α, :Δ_p_out, :γ_ns_shape, :γ_ns_on_09, :γ_ns_on_12, :η, :r,
			 :R_p, :c, :γ_s_pop, :γ_ns_pop, :s_R, :μ_R, :R_q, :γ_s_shape,
			 :γ_s_on_09, :γ_s_on_12, :σ_δ, :γ_ns_of_09_std, :βlocal]

	# Input check
	for θ_i in keys(θ_fix)
 		@assert θ_i ∈ θ_sym "Something's up: the parameter you're fixing is not in list."
 	end

	# Length of parameter vector
	N_θ = length(θ_sym)

	# Map parameter symbols to their position in parameter vector
	θ_ind = Dict(θ_sym .=> 1:N_θ)

	# Free parameters
	fix_ind  = [θ_ind[x] for x in keys(θ_fix)]
	free_ind = [θ_ind[x] for x in filter(x -> x ∉ keys(θ_fix), θ_sym)]
	fix_val  = [θ_fix[x] for x in keys(θ_fix)]

	build_θ(x::Vector{Float64}) = build_θ(x, fix_val, free_ind, fix_ind)
	obj_fun(x::Vector{Float64}) = obj(build_θ(x), distpara0, data[:on_12],
									  data[:on_09], data[:of_09]; parallel = parallel)

	x0 = deepcopy(θ_init)
	if only_likelihoods
    	return obj_fun(x0[free_ind], distpara0, data[:on_12], data[:on_09],
					   data[:of_09]; parallel = parallel, allout=true)
	end

	# bounds: x[3] < 0 || x[4] < 0 || x[12] > 1 || x[12] < 0

	# Optimize the objective function
	res = optimize(obj_fun, x0[free_ind], Optim.Options(f_calls_limit = Int(1e5),
		     	   iterations = Int(1e5), show_trace = true, store_trace = true))
	x, fval = res.minimizer, res.minimum

	CSV.write("$OUTPUT/bootstrap_estimates_$(vint).csv", Tables.table(build_θ(x)))
	return θ
end

function build_θ(θ_free_val::Vector{S}, θ_fix_val::Vector{S}, free_ind::Vector{T},
				 fix_ind::Vector{T}) where {S<:Float64, T<:Int64}
	N = length(fix_ind) + length(free_ind)
	x = zeros(N)
	x[free_ind] .= θ_free_val
	x[fix_ind]  .= θ_fix_val
	return x
end

function obj(θ::V, distpara0::V, data12::D, data09::D, bp::D;
			 parallel::Bool=true, allout = false) where {V<:Vector{Float64}, W<:Vector{Int64},
										 				 D<:Dict{Symbol,Vector{<:Number}}}
	out = try
		full_model(θ, distpara0, data12, data09, bp; parallel = parallel)
	catch err
		#print(err)
		throw(err)
		Inf
	end
	return allout ? out :  out[1]
end
