# Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>

# Set seed for testing purposes
rng = MersenneTwister(1234)
path = dirname(@__FILE__)

@show run_mode
## Input: True data, randomly generated title-level index
# ~ Contructs bootstrap dataset, runs the estimation ~
# Output: bootstrap_welfare.csv

## Estimates and welfare from true data
true_estimates = [0.32135,	5.6459,	14.855,	 1.1614,  0.6486,   1.9196,	 14.771,
                  -2.4895,	1,	    0.44004, 0.32415, 0.87235,  0.5,	 0.25921,
                  -9.1217,  80.267,	-13.647, 1.7296,  8.8188,   0.92622, 4.283,
                  4.9097,	0,      7.8609,	 7.739,	  0.011111, 6.8386,  6.5027,
                  0.028621]

## Load data
# Mat file created by Masao contains original data, random index for bootstrap.
vars = matread("$path/data/DataToRun_pop09_boot.mat")

## Data Renaming and Setup Pre-Bootstrap
γ0vec     = vcat(quantile.(Gamma(0.5, 20), 0.005:0.01:0.895), 28:2:60, 64:4:100)
δ_vec     = vcat(exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)), 3:2:20)
data12    = vars["data12nopop"]
data09    = vars["data09nopop"]
bp        = vars["bpnopop"]
bootindex = Int.(vars["bootindex"])

boot = DataFrame()
if run_mode == :EVAL
    # Mode 1 (optimize) can be run on several servers simultaneously to save time;
    # use  `unique` to remove bootstrap runs duplicated on multiple servers.
    boot = unique(CSV.read("$path/data/bootstrap_estimates.csv", DataFrame, header=false))[:,2:end]
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
    @show i
    ##  Generate bootstrap data
    bs_ind = bootindex[i,:]

    #################### 2009 ####################
    bmktsize09 = Int.(data09["mktsize"][bs_ind])
    bfirst09   = Int.(data09["first"][bs_ind])
    bcdindex09 = Int.(data09["cdindex"][bs_ind])
    mktsize_09 = Int(sum(bmktsize09))

    bdata09 = Dict{String, Any}()
    for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
              :popular, :pdif, :conditiondif, :basecond]
        bdata09[String(x)] = zeros(mktsize_09)
    end

    bend09   = Int.(cumsum(bmktsize09))
    bstart09 = vcat(1, bend09[1:end-1] .+ 1)

    for j = 1:length(bfirst09)
        rng_j_09  = bstart09[j]:bend09[j]
        rng_ij_09 = bfirst09[j]:bcdindex09[j]

        bdata09["cdid"][rng_j_09] = fill(j, length(rng_j_09), 1)
        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond]
            bdata09[String(x)][rng_j_09] .= data09[String(x)][rng_ij_09]
        end
    end
    bdata09["first"]   = bstart09'
    bdata09["cdindex"] = bend09'
    bdata09["mktsize"] = bmktsize09'
    bdata09["N"] = length(bdata09["p"])
    bdata09["M"] = length(bdata09["first"])

    #################### 2012 ####################
    bmktsize12 = data12["mktsize"][bs_ind]
    bfirst12   = data12["first"][bs_ind]
    bcdindex12 = data12["cdindex"][bs_ind]

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

    # bp
    bmktsizebp = bp["mktsize"][bs_ind]
    bfirstbp   = bp["first"][bs_ind]
    bcdindexbp = bp["cdindex"][bs_ind]
    bendbp     = cumsum(bmktsizebp)
    bstartbp   = [1; bendbp[1:end-1] .+ 1]

    bbp = Dict{String,Any}()
    bbp["cdid"] = collect(1:236)  # bp.cdid[bs_ind]
    for x in [:p, :obsweight, :numlist, :condition, :localint, :popular, :conditiondif]
        bbp[String(x)] = bp[String(x)][bs_ind]
    end
    bbp["first"]   = bstartbp'
    bbp["cdindex"] = bendbp'
    bbp["mktsize"] = bmktsizebp'
    bbp["N"]       = length(bbp["p"])
    bbp["M"]       = length(bbp["first"])

    distpara0 = vec(vars["distpara0"])

    if run_mode == :OPTIM

        x0 = true_estimates[7:20]
        objectivefun(x) = objective(x, x0, distpara0, γ0vec, δ_vec, bdata12, bdata09, bbp)
        x00 = x0
        # deleteat!(x0, [3, 7])

        res = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5), iterations = Int(1e5)))
        x, fval = res.minimizer, res.minimum

        estimates = [i; x[1:2]; x00[3]; x[3:5]; x00[7]; x[6:(length(x00)-2)]]
        CSV.write("$path/output_data/bootstrap_estimates_$(vint).csv", Tables.table(estimates))

    elseif run_mode == :EVAL

        xinitial = Vector{Float64}(boot[i,:])
        llh, newdistpara, fother, fWF = full_model(xinitial, distpara0, γ0vec,
                                        δ_vec, bdata12, bdata09, bbp; WFcal = true)
        wel09    = mean(fWF["AveWF09"])
        wel12    = mean(fWF["AveWF12"])
        weloff   = mean(fWF["WFbp"])
        result_w = [i; xinitial; newdistpara; wel09; wel12; weloff]
        CSV.write("$path/output_data/bootstrap_welfare_$(vint).csv", result_w)
    end
end
b_boot = output_statistics(; boot_out = "$path/output_data/bootstrap_welfare_$(vint).csv",
                           vint = "2022-06-26", write_out = true)[1]
make_table_results(b_boot; table_title = "estimates_20220626_mode_$(String(run_mode)).tex")
