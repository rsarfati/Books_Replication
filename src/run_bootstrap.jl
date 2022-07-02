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



		estimate_model(data = data_i, distpara0 = distpara0,
					   θ_init = θ_start[i,:], θ_fix = θ_fix, θ_lb = θ_lb,
					   θ_ub = θ_ub, eval_only = eval_only,
					   parallel = parallel, save_output = false)

		if !eval_only
			CSV.write("$path/output/data/bs_estimates_$(vint)_run=$i.csv", Tables.table([i; θ_i]))
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
"""
function index_data(d_in::Dict{Symbol, Dict{Symbol,Vector{<:Number}}}, ind::Vector{Int64})

	d = deepcopy(d_in)
	for y in keys(d)

		d[y][:mktsize] = d_in[y][:mktsize][ind]
		b_first        = d_in[y][:first][ind]
	    b_cdindex      = d_in[y][:cdindex][ind]

		mktsize = Int(    sum(d[y][:mktsize]))
		b_end   = Int.(cumsum(d[y][:mktsize]))
		b_start = [1; b_end[1:end-1] .+ 1]

		for k in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
	              :popular, :pdif, :conditiondif, :basecond]
			d[y][k] = zeros(mktsize)

			d[y][k][ind_to] .= d_in[y][ind_from]
		end

		d[y][:N] = length(d[y][:p])
		d[y][:M] = length(d[y][:first])
	end



	# bmktsize09 = Int.(data09["mktsize"][bs_ind])
    # bfirst09   = Int.(data09["first"][bs_ind])
    # bcdindex09 = Int.(data09["cdindex"][bs_ind])
    # mktsize_09 = Int(sum(bmktsize09))

    # bdata09 = Dict{String, Any}()
    # for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
    #           :popular, :pdif, :conditiondif, :basecond]
    #     bdata09[String(x)] = zeros(mktsize_09)
    # end

    # bend09   = Int.(cumsum(bmktsize09))
    # bstart09 = vcat(1, bend09[1:end-1] .+ 1)

    for j = 1:length(bfirst09)
        rng_j_09  = bstart09[j]:bend09[j]
        rng_ij_09 = bfirst09[j]:bcdindex09[j]

        bdata09["cdid"][rng_j_09] = fill(j, length(rng_j_09), 1)
        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond]
            bdata09[String(x)][rng_j_09] .= data09[String(x)][rng_ij_09]
        end
    end
    # bdata09["first"]   = bstart09'
    # bdata09["cdindex"] = bend09'
    # bdata09["mktsize"] = bmktsize09'
    # bdata09["N"] = length(bdata09["p"])
    # bdata09["M"] = length(bdata09["first"])
	return d
end


## A validation mode:
# if one wants to reproduce the true data estimates or validate the
# estimation method is the same as we used for the true data estimates, the
# following tricks can be used to reproduce the true data set.

# in mode 1:
# bootindex(1,:) = 1:236
# in mode 2, in addition to the line above:
# boot(1,:) = results2(end, 8:21)

## Run Bootstrap
for i = 1:N_bs
    ##  Generate bootstrap data
    bs_ind = bootindex[i,:]

    #################### 2009 ####################
    # bmktsize09 = Int.(data09["mktsize"][bs_ind])
    # bfirst09   = Int.(data09["first"][bs_ind])
    # bcdindex09 = Int.(data09["cdindex"][bs_ind])
    # mktsize_09 = Int(sum(bmktsize09))

    # bdata09 = Dict{String, Any}()
    # for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
    #           :popular, :pdif, :conditiondif, :basecond]
    #     bdata09[String(x)] = zeros(mktsize_09)
    # end

    # bend09   = Int.(cumsum(bmktsize09))
    # bstart09 = vcat(1, bend09[1:end-1] .+ 1)

    for j = 1:length(bfirst09)
        rng_j_09  = bstart09[j]:bend09[j]
        rng_ij_09 = bfirst09[j]:bcdindex09[j]

        bdata09["cdid"][rng_j_09] = fill(j, length(rng_j_09), 1)
        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond]
            bdata09[String(x)][rng_j_09] .= data09[String(x)][rng_ij_09]
        end
    end
    # bdata09["first"]   = bstart09'
    # bdata09["cdindex"] = bend09'
    # bdata09["mktsize"] = bmktsize09'
    # bdata09["N"] = length(bdata09["p"])
    # bdata09["M"] = length(bdata09["first"])


    #################### 2012 ####################
    # bmktsize12 = data12["mktsize"][bs_ind]
    # bfirst12   = data12["first"][bs_ind]
    # bcdindex12 = data12["cdindex"][bs_ind]
    mktsize_12 = Int(sum(bmktsize12))

    bdata12 = Dict{String, Any}()
    for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
              :popular, :pdif, :conditiondif, :basecond, :disappear]
        bdata12[String(x)] = zeros(mktsize_12)
    end

    bend12   = cumsum(bmktsize12)
    bstart12 = [1; bend12[1:end-1] .+ 1]

    for j = 1:length(bfirst12)
        j_rng_12  = Int.(bstart12[j]:bend12[j])
        ij_rng_12 = Int.(bfirst12[j]:bcdindex12[j])

        bdata12["cdid"][j_rng_12] = fill(j, length(j_rng_12), 1)

        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond, :disappear]
            bdata12[String(x)][j_rng_12] .= data12[String(x)][ij_rng_12]
        end
    end
    bdata12["first"]   = bstart12'
    bdata12["cdindex"] = bend12'
    bdata12["mktsize"] = bmktsize12'

    bdata12["N"]       = length(bdata12["p"])
    bdata12["M"]       = length(bdata12["first"])

    ##################### 2009 OFFLINE ####################
    # bmktsizebp = bp["mktsize"][bs_ind]
    # bfirstbp   = bp["first"][bs_ind]
    # bcdindexbp = bp["cdindex"][bs_ind]

	bendbp     = cumsum(bmktsizebp)
    bstartbp   = [1; bendbp[1:end-1] .+ 1]

    bbp = Dict{String,Any}()
    bbp["cdid"] = collect(1:236)  # bp.cdid[bs_ind]
    for x in [:p, :obsweight, :numlist, :condition, :localint, :popular, :conditiondif]
        bbp[String(x)] = bp[String(x)][bs_ind]
    end
    bbp["first"]   = bstartbp'
    bbp["cdindex"] = bendbp'
    # bbp["mktsize"] = bmktsizebp'
    # bbp["N"]       = length(bbp["p"])
    # bbp["M"]       = length(bbp["first"])
end
