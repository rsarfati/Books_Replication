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

# Keywords:
- `vint::String`: Specify the estimation vintage. (For storage sanity.)
- `eval_only::Bool`: Evaluate the likelihood of model for given parameter vector
	`θ_init.` Does *not* run MLE.
- `parallel::Bool`: when true, runs parts of likelihood evaluation in parallel.
- `θ_init::OrderedDict{Symbol,Float64}`: When `eval_only=true`, evaluates
	likelihood of given `θ_init`. When `eval_only=false`, runs maximum
	likelihood estimation with `θ_init` as starting parameter vector. Note that
	these values are *transformations* of the parameters in question.
- `θ_fix::Dict{Symbol,Float64}`: Maps fixed parameters to calibrated values.
- `θ_{lb,ub}::Dict{Symbol,Float64}`: Maps parameters to lower and upper bounds,
	respectively. Bounds on fixed (calibrated) parameters are ignored. (That is,
	if `θ_fix` values are in conflict with `θ_{lb,ub}`, the former takes
	precedence.)
"""
function estimate_model(; vint::String = "", eval_only = false, parallel = true,
 			θ_init::OrderedDict{Symbol,Float64} = OrderedDict(
			#=1=#	:α          	=> 14.771,  #1  α
			#=2=#	:Δ_p_out    	=> -2.4895, #2  Δ_p_out
			#=3=#	:γ_ns_shape 	=> 1.0,     #3  γ_ns_shape *** FIXED
			#=4=#	:γ_ns_on_09 	=> 0.44004, #4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
			#=5=#	:γ_ns_on_12 	=> 0.32415, #5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
			#=6=#	:η          	=> 0.87235, #6  η - 1   *** I am infering there is meant to be a -1 transformation
			#=7=#	:r          	=> 0.5,     #7  r * 10  *** FIXED
			#=8=#	:R_p        	=> 0.25921, #8  R_p
			#=9=#	:c          	=> -9.1217, #9  c * 10
			#=10=#	:γ_s_pop    	=> 80.267,  #10 γ_s_pop * 100
			#=11=#	:γ_ns_pop   	=> -13.647, #11 γ_ns_pop * 10
			#=12=#	:s_R        	=> 1.7296,  #12 s_R
			#=13=#	:μ_R        	=> 8.8188,  #13 μ_R / s_R
			#=14=#	:R_q        	=> 0.92622, #14 1 - R_q
			#=15=#	:γ_s_shape  	=> 4.283,   #15 γ_s_shape
			#=16=#	:γ_s_on_09  	=> 4.9097,  #16 γ_s_on_09
			#=17=#	:γ_s_on_12  	=> 0.0,     #17 γ_s_on_12
			#=18=#	:σ_δ        	=> 7.8609,  #18 σ_δ
			#=19=#	:γ_ns_of_09_std => 7.739,   #19 γ_ns_of_09_std
			#=20=#	:βlocal         => 0.011111), #20 βlocal ( = γ_ns_of_09_loc / γ_ns_of_09_std)
			θ_fix = Dict{Symbol,Float64}(:r => 0.5, :γ_ns_shape => 1),
			θ_lb  = Dict{Symbol,Float64}([:γ_ns_on_09, :γ_ns_on_12, :R_q] .=> 0.),
			θ_ub  = Dict{Symbol,Float64}(:R_q => 1.))

	# Load data
	@load "$INPUT/data_to_run.jld2" data
	@load "$INPUT/distpara0.jld2"   distpara0

	# In case not optimizing; just evaluate model at θ_init & return!
	if eval_only
    	return obj(vals(θ_init), distpara0, data[:on_12], data[:on_09],
				   data[:of_09]; parallel = parallel, allout = true)
	end

	# Extract parameter information
	N_θ   = length(θ_init)
	θ_sym = keys(θ_init)
	θ_val = vals(θ_init)

	# Map parameter symbols to their position in parameter vector
 	θ_ind = Dict(θ_sym .=> 1:N_θ)

	# Input check: cannot fix parameters which are not listed in θ_init
	for θ_i in [keys(θ_fix); keys(θ_lb); keys(θ_ub)]
 		@assert θ_i ∈ θ_sym "Input error: fixed/bounded parameter $θ_i unknown."
 	end

	# Extract indices of free and fixed parameters
	fix_val  = vals(θ_fix)
	fix_ind  = sort([θ_ind[x] for x in keys(θ_fix)])
	free_ind = deleteat!(collect(1:N_θ), fix_ind)

	# Function closures
	build_θ(x::Vector{Float64}) = build_θ(x, fix_val, free_ind, fix_ind)
	obj_fun(x::Vector{Float64}) = obj(build_θ(x), distpara0, data[:on_12],
									  data[:on_09], data[:of_09];
									  parallel = parallel)[1]

	# Construct vectors of lower and upper bounds for optimization routine
	lb, ub = repeat.([[-Inf], [Inf]], N_θ)
	for l in θ_lb;	lb[θ_ind[l[1]]] = l[2] end
	for u in θ_ub;	ub[θ_ind[u[1]]] = u[2] end

	# Runs optimization of objective function, then reconstitutes optimal
	# parameter vector so it again includes fixed/calibrated parameters.
	res = optimize(obj_fun, θ_val[free_ind], lb[free_ind], ub[free_ind],
				   Optim.Options(f_calls_limit=Int(1e5), iterations=Int(1e5),
		     	   show_trace = true, store_trace = true))
	θ, llh = build_θ(res.minimizer), res.minimum

	# Save output (writing CSV for legacy compatibility)
	@save     "$OUTPUT/data/estimation_results_$(vint).jld2" θ, llh
	CSV.write("$OUTPUT/data/estimation_results_$(vint).csv", Tables.table(θ))

	return θ, llh
end

"""
```
obj(θ::V, distpara0::V, data12::D, data09::D, bp::D;
	parallel::Bool=true) where {V<:Vector{Float64}, W<:Vector{Int64},
  						        D<:Dict{Symbol,Vector{<:Number}}}
```
Evaluates the model and returns the likelihood. If `allout` = false, function
will return the likelihood, alone.

# Note on Errors
If an error is thrown, first check if it's a domain error (e.g. trying to take
log/sqrt of negative number). Said error implies combination of parameters is
infeasible (assuming model), so one can conclude `θ` is "infinitely unlikely."

If *not* a domain error, something is wrong with the code! Throw it to the
console and let the user investigate the cause. :)
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
			println("To user spamming CTRL-C/D: I hear your call.")
			throw(err)
		else
			println("Mystery flavor!")
			print(err)
			throw(err)
		end
		Inf
	end
	return out
end
