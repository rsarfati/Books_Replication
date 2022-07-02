"""
```
run_bootstrap(; data::Dict = Dict(),
			    distpara0::Vector{Float64} = Vector{Float64}(),
				# Parameter specification
				θ_init::OrderedDict{Symbol,Float64} = OrderedDict(),
				θ_fix = Dict(), θ_lb  = Dict(), θ_ub  = Dict(),
				N_bs::Int64 = 200, run_mode::Symbol = :EVAL)
```
# Input: True data, randomly generated title-level indices.
# Output: bootstrap_welfare.csv

Contructs bootstrap dataset, runs the estimation.

Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>
"""
function run_bootstrap(; data::Dict = Dict(),
						 distpara0::Vector{Float64} = Vector{Float64}(),
						 # Parameter specification
						 θ_init::OrderedDict{Symbol,Float64} = OrderedDict(
						   #=1=#	:α          	=> 14.771,  # TRANSFORMED?
						   #=2=#	:Δ_p_out    	=> -2.4895,
						   #=3=#	:γ_ns_shape 	=> 1.0,
						   #=4=#	:γ_ns_on_09 	=> 0.44004, # γ_ns_on_09 * 9.5 ^ (-η) / (10 * r * γ_ns_shape)
						   #=5=#	:γ_ns_on_12 	=> 0.32415, # γ_ns_on_12 * 8.0 ^ (-η) / (10 * r * γ_ns_shape)
						   #=6=#	:η          	=> 0.87235, # η - 1  (I am infering the -1)
						   #=7=#	:r          	=> 0.5,     # r * 10
						   #=8=#	:R_p        	=> 0.25921,
						   #=9=#	:c          	=> -9.1217, # c * 10
						   #=10=#	:γ_s_pop    	=> 80.267,  # γ_s_pop * 100
						   #=11=#	:γ_ns_pop   	=> -13.647, # γ_ns_pop * 10
						   #=12=#	:s_R        	=> 1.7296,
						   #=13=#	:μ_R        	=> 8.8188,  # μ_R / s_R
						   #=14=#	:R_q        	=> 0.92622, # 1 - R_q
						   #=15=#	:γ_s_shape  	=> 4.283,
						   #=16=#	:γ_s_on_09  	=> 4.9097,
						   #=17=#	:γ_s_on_12  	=> 0.0,
						   #=18=#	:σ_δ        	=> 7.8609,
						   #=19=#	:γ_ns_of_09_std => 7.739,
						   #=20=#	:βlocal         => 0.011111), # βlocal ( = γ_ns_of_09_loc / γ_ns_of_09_std)
						 θ_fix = Dict(:r => 0.5, :γ_ns_shape => 1),
						 θ_lb  = Dict([:γ_ns_on_09, :γ_ns_on_12, :R_q] .=> 0.),
						 θ_ub  = Dict(:R_q => 1.),
						 N_bs::Int64 = 200, eval_only = false, WFcal = true,
						 seed = true, VERBOSE = true)

    # Seed for consistent output
    if seed; MersenneTwister(1234) end

	# Load bootstrap indices
    @load "$INPUT/bootstrap_indices.jld2" bootindex

	# Load data + starting values for θ and distpara0
    isempty(data)      && @load "$INPUT/data_to_run.jld2" data
	isempty(distpara0) && @load "$INPUT/distpara0.jld2"   distpara0
	θ_start = if eval_only
	    # Optimization can be run on several servers simultaneously to save time;
	    # we use `unique` to remove bootstrap runs duplicated on multiple servers.
	    unique(readdlm("$INPUT/bootstrap_estimates.csv", ','), dims=1)[:,2:end]
	else
		# If running optimizer, initial values are last estimation's results
		repeat(vals(θ_init)', N_bs)
	end

	# Run bootstrap! TODO: parallelization option
    for i=1:N_bs
		println(VERBOSE, "Bootstrap iteration: $i")
		bs_ind = bootindex[i,:]

		# TODO: indexing!
		data_i = index_data(data, bs_ind)


		estimate_model(data = data_i, distpara0 = distpara0,
					   θ_init = θ_start[i,:], θ_fix = θ_fix, θ_lb = θ_lb,
					   θ_ub = θ_ub, eval_only = eval_only,
					   parallel = parallel, save_output = false)

		if !eval_only
			CSV.write("$OUTPUT/data/bs_estimates_$(vint)_run=$i.csv", Tables.table([i; θ_i]))
		end
	end

	# Output results from estimation
	b_boot = output_statistics(; boot_out = "$OUTPUT/data/bootstrap_welfare_$(vint).csv",
	                             vint = "2022-06-26", write_out = true)[1]
	make_table_results(b_boot;
		table_title = "$OUTPUT/tables/estimates_$(vint)_eval_only=$(eval_only).tex")

	if WFcal
		wel09    = mean(fWF["AveWF09"])
		wel12    = mean(fWF["AveWF12"])
		weloff   = mean(fWF["WFbp"])
		result_w = [i; xinitial; newdistpara; wel09; wel12; weloff]
		CSV.write("$OUTPUT/data/bs_welfare_$(vint).csv", result_w)
	end
	return b_boot
end

"""
```
index_data(d_in::Dict{Symbol, Dict{Symbol,Vector{<:Number}}}, ind::Vector{Int64})
```
Index into data dictionary for a given set of bootstrapped indices `ind`.
Returns bootstrapped data in a dictionary of the same necessary input format.
"""
function index_data(d_in::Dict{Symbol,Dict{Symbol,Vector{<:Number}}}, ind::Vector{Int64})

	d = deepcopy(d_in)

	obs_list = Dict{Symbol,Vector{Symbol}}(:of_09 =>
		 								   [:p, :obsweight, :numlist, :condition,
										   	:localint, :popular, :conditiondif])
	obs_list[:on_09] = [obs_list[:of_09]; :pdif; :basecond]
	obs_list[:on_12] = [obs_list[:on_09]; :disappear]

	for y in keys(d)

		d[y][:mktsize] = d_in[y][:mktsize][ind]
		b_first        = d_in[y][:first][ind]
	    b_cdindex      = d_in[y][:cdindex][ind]

		mktsize =    sum(d[y][:mktsize])
		b_end   = cumsum(d[y][:mktsize])
		b_start = [1; b_end[1:end-1] .+ 1]

		for k in [obs_list[y]; :cdid]
			d[y][k] = zeros(mktsize)
		end

		if y ∈ [:on_09, :on_12]
			for j=1:length(ind)
				ind_to   = b_start[j]:b_end[j]
				ind_from = b_first[j]:b_cdindex[j]

				d[y][:cdid][ind_to] = fill(j, length(ind_to), 1)

				for k in obs_list[y]
					d[y][k][ind_to] .= d_in[y][ind_from]
				end
			end
		else
			d[y][:cdid] = collect(1:236) # TODO: this commented: bp.cdid[bs_ind]
			for k in obs_list[y]
				d[y][k] = d_in[y][k][ind]
			end
	    end

		d[y][:first]   = b_start
		d[y][:cdindex] = b_end

		d[y][:N] = length(d[y][:p])
		d[y][:M] = length(d[y][:first])
	end
	return d
end
