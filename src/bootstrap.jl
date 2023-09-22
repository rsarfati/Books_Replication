"""
```
run_bootstrap(; data::Dict = Dict(),
			    distpara0::Vector{Float64} = Vector{Float64}(),
				# Parameter specification
				θ_init::OrderedDict{Symbol,Float64} = OrderedDict(),
				θ_fix = Dict(), θ_lb  = Dict(), θ_ub  = Dict(),
				N_bs::Int64 = 200, run_mode::Symbol = :EVAL)
```
# Input: True data and randomly-generated title-level indices.
# Output: bootstrap_welfare.csv

Contructs bootstrap dataset, runs the estimation.

Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>
"""
function run_bootstrap(; data::Dict = Dict(),
						 distpara0::Vector{Float64} = Vector{Float64}(),
						 # Parameter specification
						 θ_init::OrderedDict{Symbol,T} = OrderedDict(),
						 θ_fix::Dict{Symbol,T} = Dict(:r => 0.5, :γ_ns_shape => 1.0),
  						 θ_lb::Dict{Symbol,T}  = Dict([:γ_ns_on_09, :γ_ns_on_12, :R_q, :R_p] .=> 0.),
  						 θ_ub::Dict{Symbol,T}  = Dict([:R_q, :R_p] .=> 1.),
						 # Specification: options are :standard, :condition, :cond_list
						 spec::Symbol = :standard,
						 # Bootstrap options
						 N_bs::Int64 = 100, seed::Bool = true,
						 bs_inds::UnitRange{Int64} = 1:N_bs,
						 vint::String  = "latest", write_output::Bool = true,
						 eval_only::Bool = false, WFcal::Bool = false,
						 parallel::Bool  = true, read_draws::String = "") where T<:Float64

    # Seed for consistent output
    if seed; MersenneTwister(1234) end

	# Overwrite if specific range is provided
	# (for parallelizing over main launch scripts)
	N_bs = length(bs_inds)

	# Load bootstrap indices
    @load "$INPUT/bootstrap_indices.jld2" bootindex

	# Load data + starting values for θ and distpara0
    isempty(data)      && @load "$INPUT/data_to_run.jld2" data
	isempty(distpara0) && @load "$INPUT/distpara0.jld2"   distpara0

	θ_init = if read_draws != ""
	    # Optimization can be run on several servers simultaneously to save time;
	    # we use `unique` to remove bootstrap runs duplicated on multiple servers.
		# NOT CURRENTLY WORKING.
		unique(readdlm(read_draws, ','), dims=1)[:,2:end]
		# To read in old draws:
	    #unique(readdlm("$INPUT/bootstrap_estimates.csv", ','), dims=1)[:,2:end]
		OrderedDict(keys(θ_init) .=> θ_start[i,:])
	end

	# Run bootstrap!
    for i=bs_inds
		println(VERBOSE, "Bootstrap iteration: $i")
		llh_i, θ_i, distpara_i = estimate_model(data = index_data(data, bootindex[i,:]),
					   		    	θ_init = θ_init, distpara0 = distpara0,
					    	    	θ_fix = θ_fix, θ_lb = θ_lb, θ_ub = θ_ub,
									spec = spec, eval_only = eval_only, WFcal = false,
					    	    	parallel = parallel, write_output = false,
									bootstrap = true)
		#if write_output
			@save     "$OUTPUT/estimation_results_$(string(spec))_$(vint)_run=$i.jld2"  θ llh distpara
			CSV.write("$OUTPUT/bs_llh_theta_$(string(spec))_$(vint)_run=$i.csv", Tables.table([llh_i; θ_i]))
			CSV.write("$OUTPUT/bs_distpara_$(string(spec))_$(vint)_run=$i.csv",  Tables.table(distpara_i))
		#end
		println("=======================================")
		println("Completed Bootstrap iteration = $(i)!")
		println("=======================================")
	end

	# if WFcal
	# 	wel09    = mean(fWF["AveWF09"])
	# 	wel12    = mean(fWF["AveWF12"])
	# 	weloff   = mean(fWF["WFbp"])
	# 	result_w = [i; xinitial; newdistpara; wel09; wel12; weloff]
	# 	CSV.write("$OUTPUT/bs_welfare_$(vint).csv", result_w)
	# end
end

"""
```
index_data(d_in::Dict{Symbol, Dict{Symbol,Vector{<:Number}}}, ind::Vector{Int64})
```
Index into data dictionary for a given set of bootstrapped indices `ind`.
Returns bootstrapped data in a dictionary of the same necessary input format.
"""
function index_data(d_in::Dict{Symbol,Dict{Symbol,Vector{<:Number}}}, ind::Vector{Int64})

	# Initialize an empty dictionary
	d = Dict(keys(d_in) .=> [Dict{Symbol,Vector{<:Number}}(),
							 Dict{Symbol,Vector{<:Number}}(),
							 Dict{Symbol,Vector{<:Number}}()])

	# Each online/offline x year combination doesn't have same set of obs.
	obs_list = Dict(:of_09 => [:p, :obs_w, :numlist, :condition,
							   :localint, :popular, :cond_dif])
	obs_list[:on_09] = [obs_list[:of_09]; :pdif; :basecond]
	obs_list[:on_12] = [obs_list[:on_09]; :disappear]

	for y in keys(d_in)
		d[y][:mktsize] = d_in[y][:mktsize][ind]
		b_first        = d_in[y][:d_first][ind]
	    b_cdindex      = d_in[y][:cdindex][ind]

		mktsize =    sum(d[y][:mktsize])
		b_end   = cumsum(d[y][:mktsize])
		b_start = [1; b_end[1:end-1] .+ 1]

		for k in [obs_list[y]; :cdid]
			d[y][k] = zeros(eltype(d_in[y][k]), mktsize)
		end

		if y ∈ [:on_09, :on_12]
			for j=1:length(ind)
				ind_to   = b_start[j]:b_end[j]
				ind_from = b_first[j]:b_cdindex[j]

				d[y][:cdid][ind_to] = fill(j, length(ind_to))
				for k in obs_list[y]
					d[y][k][ind_to] .= d_in[y][k][ind_from]
				end
			end
		else
			d[y][:cdid] = collect(1:236) # TODO: was commented bp.cdid[bs_ind]
			for k in obs_list[y]
				d[y][k] = d_in[y][k][ind]
			end
	    end
		d[y][:d_first] = deepcopy(b_start)
		d[y][:cdindex] = deepcopy(b_end)
	end
	return d
end
