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
 					  	  θ_init::OrderedDict{Symbol,T} = OrderedDict(),
						  θ_fix::Dict{Symbol,T} = Dict(:r => 0.5, :γ_ns_shape => 1.0),
						  θ_lb::Dict{Symbol,T}  = Dict([:γ_ns_on_09, :γ_ns_on_12, :R_q, :R_p] .=> 0.),
						  θ_ub::Dict{Symbol,T}  = Dict([:R_q, :R_p] .=> 1.),
						  # Specification: options are :standard, :condition, :cond_list
						  spec::Symbol = :standard,
						  # Options
						  vint::String    = "latest", write_output::Bool = true,
						  eval_only::Bool = false, WFcal::Bool = false,
						  parallel::Bool  = true, bootstrap::Bool = false) where T<:Float64

	# Load data if not provided at function call
	isempty(data)      && @load "$INPUT/data_to_run.jld2" data
	isempty(distpara0) && @load "$INPUT/distpara0.jld2"   distpara0
	if isempty(θ_init)
		θ_init = OrderedDict(
		#=1=#	:α          	=> 14.7544,	#1  α 				[15.77 (1.28)]
		#=2=#	:Δ_p_out    	=> -2.5128,	#2  Δ_p_out 		[0.16 (0.02)]
		#=3=#	:γ_ns_shape 	=> 1.0,		#3  γ_ns_shape *** FIXED
		#=4=#	:γ_ns_on_09 	=> 0.4247,	#4  γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
		#=5=#	:γ_ns_on_12 	=> 0.3275,	#5  γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
		#=6=#	:η          	=> 0.8302,	#6  η - 1 			[1.87 (0.08)]
		#=7=#	:r          	=> 0.5,     #7  r * 10     *** FIXED
		#=8=#	:R_p        	=> 0.2738,	#8  R_p 			[0.26 (0.02)]
		#=9=#	:c          	=> -9.0273,	#9  c * 10 			[-0.91 (0.07)]*10
		#=10=#	:γ_s_pop    	=> 78.7921,	#10 γ_s_pop * 100 	[0.80 (0.09)]
		#=11=#	:γ_ns_pop   	=> -14.592,	#11 γ_ns_pop * 10 	[-1.36 (0.19)]
		#=12=#	:s_R        	=> 1.7247,	#12 s_R 			[1.73 (0.09)]
		#=13=#	:μ_R        	=> 8.8787,  #13 μ_R / s_R		[15.25 (0.80) / 1.73 (0.09)]
		#=14=#	:R_q        	=> 0.9253	#14 1 - R_q       1-[0.07 (0.01)]
		# Distributional parameters (pinned down within model -- do not optimize over!)
		##=15=#	:γ_s_shape  	=> 0.3598,	#15 γ_s_shape 		[0.32 (0.03)]
		##=16=#	:γ_s_on_09  	=> 5.634,	#16 γ_s_on_09 		[5.65 (1.41)]
		##=17=#	:γ_s_on_12  	=> 13.2618,	#17 γ_s_on_12 		[14.86 (3.10)]
		##=18=#	:σ_δ        	=> 1.16,	#18 σ_δ 			[1.16 (0.05)]
		##=19=#	:γ_ns_of_09_std => 0.65,	#19 γ_ns_of_09_std 	[0.65 (0.17)]
		##=20=#	:βlocal			=> 1.923	#20 βlocal (= γ_ns_of_09_loc / γ_ns_of_09_std) [1.25 (0.37) / 0.65 (0.17)]
		)
	end
	# Create indicator for having the lowest price / "being listed first"
	for d in [data[:on_09], data[:on_12], data[:of_09]]
		d[:has_min_p] = zeros(length(d[:cdid]))
		# Iterate over 236 titles in each sample
		for t in unique(d[:cdid])
			# Find the minimum price among all listings of title t
			p_min_t = minimum(d[:p][d[:cdid] .== t])
			# Count how many listings offered title t at said minimum price
			N_min_p = sum(d[:p][d[:cdid] .== t] .== p_min_t)
			# For those listings with minimum price, assign indicator
			d[:has_min_p][(d[:cdid] .== t) .& (d[:p] .== p_min_t)] .= 1 / N_min_p
		end
	end

	@everywhere data      = @eval $data
	@everywhere distpara0 = @eval $distpara0
	@everywhere d_on_09   = data[:on_09]
	@everywhere d_on_12   = data[:on_12]
	@everywhere bp        = data[:of_09]

	# Load specification
	@assert spec ∈ [:standard, :condition, :cond_list]
	if spec != :standard
		!haskey(θ_init, :η_c) && (θ_init[:η_c] = 0.16) #=21=#
		!haskey(θ_init, :α_c) && (θ_init[:α_c] = 5.00) #=22=#
	end
	if spec == :cond_list
		!haskey(θ_init, :min_p) && (θ_init[:min_p] = 0.20) #=23=#
	end

	# Simply evaluate likelihood at θ_init & return
	if eval_only
    	out = obj(vals(θ_init), distpara0, data[:on_12], data[:on_09], data[:of_09];
				  WFcal = WFcal, parallel = parallel, spec = spec)

		println("Evaluated model at θ = ", vals(θ_init), "\nLLH:  $(out[1])")

		if bootstrap; return out[1:3] end

		if write_output
			@save "$OUTPUT/evaluated_results_$(string(spec))_$(vint).jld2" out
		end
		return out
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

		out = obj(θ_full(x), distpara0, data[:on_12], data[:on_09], data[:of_09];
				  parallel = parallel, spec = spec)

		println("LLH: $(out[1])")
		return out[1]
	end

	# Construct vectors of lower and upper bounds for optimization routine
	lb, ub = repeat.([[-Inf], [Inf]], N_θ)
	for l in θ_lb;	lb[θ_ind[l[1]]] = l[2] end
	for u in θ_ub;	ub[θ_ind[u[1]]] = u[2] end

	# Optimize objective function, then reconstitute optimal parameter
	# vector to again include fixed/calibrated parameters.
	res = optimize(obj_fun, lb[free_ind], ub[free_ind],
				   θ_val[free_ind], Fminbox(),#NelderMead(),
				   Optim.Options(f_tol = 1e-2, f_calls_limit = Int(1e4),
		     	   show_trace = true, store_trace = true))
	θ, llh = θ_full(res.minimizer), res.minimum

	out = full_model(θ, distpara0, data[:on_12], data[:on_09], data[:of_09];
			   		 spec = spec, WFcal = WFcal, parallel = parallel)
	distpara = out[3]

	# Save output (writing to CSV for legacy compatibility)
	if write_output
		@save     "$OUTPUT/estimation_results_$(string(spec))_$(vint).jld2"  θ llh distpara
		CSV.write("$OUTPUT/estimation_theta_$(string(spec))_$(vint).csv",    Tables.table(θ))
		CSV.write("$OUTPUT/estimation_distpara_$(string(spec))_$(vint).csv", Tables.table(distpara))
	end
	return llh, θ, distpara
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
function obj(θ::V, distpara0::V, data12::D, data09::D, bp::D; WFcal = false,
			 parallel = true, spec = :standard) where {V<:Vector{Float64}, W<:Vector{Int64},
									                   D<:Dict{Symbol,Vector{<:Number}}}
	out = try
		full_model(θ, distpara0, data12, data09, bp;
				   spec = spec, WFcal = WFcal, parallel = parallel)
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
