"""
```
estimate_model(; vint::String = "", eval_only = false, parallel = true,
				 θ_init::OrderedDict{Symbol,Float64} = OrderedDict()
			   	 θ_fix = Dict{Symbol,Float64}(),
			   	 θ_lb  = Dict{Symbol,Float64}(),
			   	 θ_ub  = Dict{Symbol,Float64}())
```
Estimates model given full dataset, starting from initial parameter values, as
specified by user.

# Data specification
- `data::Dict`: Provide data sample one would like to estimate model with
	respect to. When left blank, loads full dataset!
- `distpara0::Vector{Float64}`: Provide input `distpara0`, otherwise read in
	default vector.

# Parameter specification
- `θ_init::OrderedDict{Symbol,Float64}`: When `eval_only=true`, evaluates
	likelihood of given `θ_init`. When `eval_only=false`, runs maximum
	likelihood estimation with `θ_init` as starting parameter vector. Note that
	these values are *transformations* of the parameters in question.
- `θ_fix::Dict{Symbol,Float64}`: Maps fixed parameters to calibrated values.
- `θ_{lb,ub}::Dict{Symbol,Float64}`: Maps parameters to lower and upper bounds,
	respectively. Bounds on fixed (calibrated) parameters are ignored. (That is,
	if `θ_fix` values are in conflict with `θ_{lb,ub}`, the former takes
	precedence.)

# Options
- `vint::String`: Specify the estimation vintage. (For storage sanity.)
- `eval_only::Bool`: Evaluate the likelihood of model for given parameter vector
	`θ_init.` Does *not* run MLE.
- `write_output::Bool`: Do you want to save output of this estimation to file?
- `parallel::Bool`: when true, runs parts of likelihood evaluation in parallel.
"""
function estimate_model(; # Data specification
						  data::Dict = Dict(),
						  distpara0::Vector{T} = Vector{Float64}(),
						  # Parameter specification
 					  	  θ_init::OrderedDict{Symbol,T} = OrderedDict(
							#=1=#	:α          	=> 16.1319,#14.771,  #1  α
							#=2=#	:Δ_p_out    	=> 0.083,#-2.4895, #2  Δ_p_out
							#=3=#	:γ_ns_shape 	=> 1.0,     #3  γ_ns_shape *** FIXED
							#=4=#	:γ_ns_on_09 	=> 0.40863, #4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
							#=5=#	:γ_ns_on_12 	=> 0.2902, #5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
							#=6=#	:η          	=> 0.9405,#0.87235, #6  η - 1   *** I am infering there is meant to be a -1 transformation
							#=7=#	:r          	=> 0.5,     #7  r * 10  *** FIXED
							#=8=#	:R_p        	=> 0.926,#0.25921, #8  R_p
							#=9=#	:c          	=> -4.84,#-9.1217, #9  c * 10
							#=10=#	:γ_s_pop    	=> 63.948,#80.267,  #10 γ_s_pop * 100
							#=11=#	:γ_ns_pop   	=> -8.224, #11 γ_ns_pop * 10
							#=12=#	:s_R        	=> 1.833,#1.7296,  #12 s_R
							#=13=#	:μ_R        	=> 8.8417,  #13 μ_R / s_R
							#=14=#	:R_q        	=> 0.9041, #14 1 - R_q
							#=15=#	:γ_s_shape  	=> 0.367,#4.283,   #15 γ_s_shape
							#=16=#	:γ_s_on_09  	=> 5.545,#4.9097,  #16 γ_s_on_09
							#=17=#	:γ_s_on_12  	=> 14.487,#0.0,     #17 γ_s_on_12
							#=18=#	:σ_δ        	=> 7.9289,  #18 σ_δ
							#=19=#	:γ_ns_of_09_std => 6.918,   #19 γ_ns_of_09_std
							#=20=#	:βlocal			=> 0.0156, #20 βlocal (= γ_ns_of_09_loc / γ_ns_of_09_std)
							#=21=#	:η_c		    => 0.136,
							#=22=#	:α_c			=> 0.146),
						  θ_fix::Dict{Symbol,T} = Dict(:r => 0.5, :γ_ns_shape => 1.0),
						  θ_lb::Dict{Symbol,T}  = Dict([:Δ_p_out, :γ_ns_on_09, :γ_ns_on_12, :R_q] .=> 0.),
						  θ_ub::Dict{Symbol,T}  = Dict(:R_q => 1.),
						  # Options
						  vint::String    = "", write_output::Bool = true,
						  eval_only::Bool = false,
						  parallel::Bool  = true) where T<:Float64

	# Load data if not provided at function call
	isempty(data)      && @load "$INPUT/data_to_run.jld2" data
	isempty(distpara0) && @load "$INPUT/distpara0.jld2"   distpara0

	# Simply evaluate likelihood at θ_init & return
	if eval_only
    	return obj(vals(θ_init), distpara0, data[:on_12], data[:on_09], data[:of_09];
				   parallel = parallel)
	end

	# Extract information on parameters
	N_θ   = length(θ_init)
	θ_sym = keys(θ_init)
	θ_val = vals(θ_init)

	# Map parameter symbols to their position in parameter vector
 	θ_ind = Dict(θ_sym .=> 1:N_θ)

	# Input check: cannot fix parameters which are not listed in θ_init
	for θ_i in union(keys(θ_fix), keys(θ_lb), keys(θ_ub))
 		@assert θ_i ∈ θ_sym "Input error: fixed/bounded parameter $θ_i unknown."
 	end

	# Extract indices of free and fixed parameters
	fix_val  = vals(θ_fix)::Vector{Float64}
	fix_ind  = sort([θ_ind[x] for x in keys(θ_fix)])
	free_ind = deleteat!(collect(1:N_θ), fix_ind)

	# Build function closures for optimizer
	θ_full(x::Vector{Float64}) = build_θ(x, fix_val, free_ind, fix_ind)
	function obj_fun(x::Vector{Float64})
		θ = θ_full(x)
		# Parameter is "infinitely unlikely" if out of bounds
		for l in θ_lb; if θ[θ_ind[l[1]]] < l[2]; return Inf end end
		for u in θ_ub; if θ[θ_ind[u[1]]] > u[2]; return Inf end end
		println("Parameters in-bounds, θ: $θ")

		out = obj(θ_full(x), distpara0, data[:on_12],
			  data[:on_09], data[:of_09]; parallel = parallel)

    	println("LLH: $(out[1])")
		return out[1]
	end

	# Construct vectors of lower and upper bounds for optimization routine
	#lb, ub = repeat.([[-Inf], [Inf]], N_θ)
	#for l in θ_lb;	lb[θ_ind[l[1]]] = l[2] end
	#for u in θ_ub;	ub[θ_ind[u[1]]] = u[2] end

	# Optimize objective function, then reconstitute optimal parameter
	# vector to again include fixed/calibrated parameters.
	res = optimize(obj_fun, #lb[free_ind], ub[free_ind],
				   θ_val[free_ind], Newton(),#Fminbox(),
				   Optim.Options(f_calls_limit = Int(1e5), iterations = Int(1e5),
		     	   show_trace = true, store_trace = true))
	θ, llh = θ_full(res.minimizer), res.minimum

	# Save output (writing to CSV for legacy compatibility)
	if write_output
		@save     "$OUTPUT/estimation_results_$(vint).jld2" θ llh
		CSV.write("$OUTPUT/estimation_results_$(vint).csv", Tables.table(θ))
	end
	return θ, llh
end

"""
```
obj(θ::V, distpara0::V, data12::D, data09::D, bp::D;
	parallel::Bool=true) where {V<:Vector{Float64}, W<:Vector{Int64},
  						        D<:Dict{Symbol,Vector{<:Number}}}
```
Evaluates model and returns likelihood. Catches errors implying combination
of parameters is infeasible under model (e.g. domain error arising from taking
log/sqrt of negative number), and returns "infinitely unlikely."
"""
function obj(θ::V, distpara0::V, data12::D, data09::D, bp::D;
			 parallel = true) where {V<:Vector{Float64}, W<:Vector{Int64},
									 D<:Dict{Symbol,Vector{<:Number}}}
	out = try
		full_model(θ, distpara0, data12, data09, bp; parallel = parallel)
	catch err
		if typeof(err)<:DomainError
			println("Domain error in optimization!")
			print(err)
		elseif typeof(err)<:InterruptException
			println("\n... To the user spamming CTRL-C/D: I hear your call.")
			throw(err)
		else
			println("Mystery flavor!")
			print(err)
			throw(err)
		end
		Inf, Vector{Float64}(), Dict{String,Any}(), Dict{String,Any}(), Inf, Inf
	end
	return out
end
