# Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>
using CSV, DataFrames, Dates, Distributions, MAT, Optim, SparseArrays, Statistics

## Input: True data, randomly generated title-level index
# ~ Contructs bootstrap dataset, runs the estimation ~
# Output: bootstrap_welfare.csv

## Important notes

# VARIABLE NAMES
# - fullmodelllhWFAug22newtest2015 -> full_model
# - fullmodelllhWFAug22newtest2015_all -> full_model_all
# - current_mode -> mode

# RELATED FILES
# standard_errors.m (by Masao) processed <bootstrap_welfare.csv> to generate
# bootstrap numbers in Summary201609.xls. Said code is now appended to the end
# of this file, so there is no longer a need to run standard_errors.m

## Two modes
# This file has two modes. Mode 2 can be run after Mode 1.

# In the first mode, it cold starts to estimates parameters in the first
# input of the likelihood function full_model_all with
# bootstrap sample, and save the results in bootstrap_estimates.csv. These
# parameters are/transform to alpha, beta, gammaishape, gammaimean09,
# gammaimean12, eta, r, olp, c, lamda1, lamda2, betacond, betapop, olm,
# oltheta, naturaldisappear. This mode is very computationally intensive
# (involves optimization) and has to be run on servers and still takes
# weeks.

# In the second mode, it figures out the six parameters and welfare numbers
# that full_model_all can compute internally, reading
# the estimates from the first mode as input. This mode is quite fast
# without the welfare part, and can be run on a desktop in a couple hours.
# But it will take a day on servers if welfare output (WF) from
# full_model_all is included.

# As of Dec 2017, Mode 1 is done with the issue of duplicated book
# index in bootstrap sample, Mode 2 is done with the issue
# corrected, albeit taking Mode's 1's results as input.

## Estimates and welfare from true data
true_estimates = [0.32135,	5.6459,	14.855,	 1.1614,  0.6486,   1.9196,	 14.771,
                  -2.4895,	1,	    0.44004, 0.32415, 0.87235,  0.5,	 0.25921,
                  -9.1217,  80.267,	-13.647, 1.7296,  8.8188,   0.92622, 4.283,
                  4.9097,	0,      7.8609,	 7.739,	  0.011111, 6.8386,  6.5027,
                  0.028621]
# (These are the estimates and welfare from true data. It is copied from last row
# of results2.csv created by Masao. Except for the 5th and 6th elements,
# E[gamma_i] 2009 offline and betalocal, which according to Hongkai's validation
# mode, is different from Masao's number.)

# TODO: Specify script parameters
vint      = "2022-06-04"
run_mode  = 2   # Choose between available modes: 1 or 2 (see description above)
N_workers = 72  # No. workers to request from cluster
N_bs      = 200 # No. bootstrap iterations

## Parallel computing setup
# The following 3 lines are for reduced-form server.
# Can change resource parameters (walltime and mem) if necessary
# dcluster = parcluster
# dcluster.ResourceTemplate = '-l nodes=^N^, software = MATLAB_Distrib_Comp_Engine+^N^,
# walltime = 80:00:00, mem = 64gb'
# dcluster.saveProfile

# This line is common for all servers or desktop. Choose number of workers! accordingly.
#parpool(N_workers)

## Load data
# Mat file created by Masao contains original data, random index for bootstrap.
vars = matread("data/DataToRun_pop09_boot.mat")

## Data Renaming and Setup Pre-Bootstrap
γ0vec = vcat(quantile.(Gamma(0.5, 20), 0.005:0.01:0.895), 28:2:60, 64:4:100)
δ_vec  = vcat(exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)), 3:2:20)
data12    = vars["data12nopop"]
data09    = vars["data09nopop"]
bp        = vars["bpnopop"]
bootindex = Int.(vars["bootindex"])

bmktsize09 = bootindex
bfirst09   = bootindex
bcdindex09 = bootindex

bmktsize12 = bootindex
bfirst12   = bootindex
bcdindex12 = bootindex

bmktsizebp = bootindex
bfirstbp   = bootindex
bcdindexbp = bootindex

# Meta Programming
# f1 = eval(Meta.parse("@formula(xt_lag0 ~ " *
#                      (*).(["$(s_p(p))xt1_lag$(p)" for p=0:lags-1]...) *
#                      (lags>1 ? (*).([" + Dzt_lag$(p)" for p=1:lags-1]...) : "") *
#                      " + (Dzt_lag0  ~ " *
#                      (*).(["$(s_p(p))zt1_lag$(p)" for p=0:lags-1]...) * "))"))

if run_mode == 2
    # Mode 1 can be run on several servers simultaneously to save time;
    # use  `unique` to remove bootstrap runs duplicated on multiple servers.
    boot = unique(CSV.read("data/bootstrap_estimates.csv",
                           DataFrame, header=false))[:,2:end]
end

# A validation mode:
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
    bmktsize09[i,:] = Int.(data09["mktsize"][bs_ind])
    bfirst09[i,:]   = Int.(data09["first"][bs_ind])
    bcdindex09[i,:] = Int.(data09["cdindex"][bs_ind])
    mktsize_09      = Int(sum(bmktsize09[i,:]))

    bdata09 = Dict{String, Any}()
    for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
              :popular, :pdif, :conditiondif, :basecond]
        bdata09[String(x)] = zeros(mktsize_09)
    end

    bend09 = Int.(cumsum(bmktsize09[i,:]))
    bstart09 = vcat(1, bend09[1:end-1] .+ 1)

    for j = 1:length(bfirst09[i,:])
        rng_j_09  = bstart09[j]:bend09[j]
        rng_ij_09 = Int.(bfirst09[i,j]:bcdindex09[i,j])

        bdata09["cdid"][rng_j_09] = fill(j, length(rng_j_09), 1)
        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond]
            bdata09[String(x)][rng_j_09] .= data09[String(x)][rng_ij_09]#
        end
    end
    bdata09["first"]   = bstart09'
    bdata09["cdindex"] = bend09'
    bdata09["mktsize"] = bmktsize09[i,:]'
    bdata09["N"] = length(bdata09["p"])
    bdata09["M"] = length(bdata09["first"])

    #################### 2012 ####################
    bmktsize12[i,:] = data12["mktsize"][bs_ind]
    bfirst12[i,:]   = data12["first"][bs_ind]
    bcdindex12[i,:] = data12["cdindex"][bs_ind]

    mktsize_12 = Int(sum(bmktsize12[i,:]))
    bdata12 = Dict{String, Any}()
    for x in [:p, :cdid, :obsweight, :numlist, :condition, :localint,
              :popular, :pdif, :conditiondif, :basecond, :disappear]
        bdata12[String(x)] = zeros(mktsize_12)
    end

    bend12   = cumsum(bmktsize12[i,:])
    bstart12 = vcat(1, bend12[1:end-1] .+ 1)

    for j = 1:length(bfirst12[i,:])
        j_rng_12  = Int.(bstart12[j]:bend12[j])
        ij_rng_12 = Int.(bfirst12[i,j]:bcdindex12[i,j])

        bdata12["cdid"][j_rng_12] = fill(j, length(j_rng_12), 1)

        for x in [:p, :obsweight, :numlist, :condition, :localint,
                  :popular, :pdif, :conditiondif, :basecond, :disappear]
            bdata12[String(x)][j_rng_12] .= data12[String(x)][ij_rng_12]
        end
    end
    bdata12["first"]   = bstart12'
    bdata12["cdindex"] = bend12'
    bdata12["mktsize"] = bmktsize12[i,:]'
    bdata12["N"]       = length(bdata12["p"])
    bdata12["M"]       = length(bdata12["first"])

    # bp
    bmktsizebp[i,:] = bp["mktsize"][bs_ind]
    bfirstbp[i,:]   = bp["first"][bs_ind]
    bcdindexbp[i,:] = bp["cdindex"][bs_ind]
    bendbp          = cumsum(bmktsizebp[i,:])
    bstartbp        = vcat(1, bendbp[1:end-1] .+ 1)

    bbp = Dict{String,Any}()
    bbp["cdid"] = collect(1:236)  # bp.cdid[bs_ind]
    for x in [:p, :obsweight, :numlist, :condition, :localint, :popular, :conditiondif]
        bbp[String(x)] = bp[String(x)][bs_ind]
    end
    bbp["first"]   = bstartbp'
    bbp["cdindex"] = bendbp'
    bbp["mktsize"] = bmktsizebp[i,:]'
    bbp["N"]       = length(bbp["p"])
    bbp["M"]       = length(bbp["first"])

    # Optimization
    if run_mode == 1

        x0 = true_estimates[7:20]
        objectivefun(x) = objective(x, x0, vec(vars["distpara0"]), γ0vec, δ_vec, bdata12, bdata09, bbp)
        x00 = x0
        # TODO: investigate why these columns are being deleted
        #deleteat!(x0, [3, 7])

        res = optimize(objectivefun, x0, Optim.Options(f_calls_limit = Int(1e5), iterations = Int(1e5)))
        x, fval = res.minimizer, res.minimum

        estimates = [i, x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
        CSV.write("bootstrap_estimates_$(vint).csv", Tables.table(estimates))

    elseif run_mode == 2

        xinitial = boot[i,:]
        llh, newdistpara, fother, fWF = full_model(xinitial, distpara0, γ0vec, δ_vec, bdata12, bdata09, bbp; WFcal = true)

        wel09    = mean(fWF["AveWF09"])
        wel12    = mean(fWF["AveWF12"])
        weloff   = mean(fWF["WFbp"])
        result_w = [i, xinitial, newdistpara, wel09, wel12, weloff]

        CSV.write("bootstrap_welfare.csv", result_w)
    end
end

## Lines adapted from Masao's standard_errors.m, can be run after Mode 2.

# Bootstrap β_σ statistics
boot = CSV.read("data/bootstrap_welfare.csv", DataFrame, header=false)
boot = unique(boot)
boot = boot[:, 2:15]

distpara = CSV.read("data/bootstrap_welfare.csv", DataFrame, header=false)
distpara = unique(distpara)
distpara = distpara[:, 16:21]

b_boot = zeros(size(boot,1), 25)
for i = 1:size(boot,1)
    xx = Vector(boot[i,:])
    est = vcat(Vector(distpara[i,:]), xx[1], xx[2]/(1+xx[1]), xx[3],
               xx[4]*10*xx[7]/10/9.5^(-xx[6]-1), xx[5]*10*xx[7]/10/8^(-xx[6]-1),
               xx[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1], 0, 0, xx[12:13], xx[14])
    b_boot[i,1] = est[1]
    b_boot[i,2] = est[2]
    b_boot[i,3] = est[3]
    b_boot[i,4] = est[7]
    b_boot[i,5] = est[8]
    b_boot[i,6] = est[13]
    b_boot[i,7] = est[15]
    b_boot[i,8] = est[4]
    b_boot[i,9] = est[12]
    b_boot[i,10] = est[9]
    b_boot[i,11] = est[10] .* est[9]
    b_boot[i,12] = (est[10] .* est[9] ./ 10) .^ (est[12] .+ 1)
    b_boot[i,13] = est[11] .* est[9]
    b_boot[i,14] = (est[11] .* est[9] ./ 10) .^ (est[12] .+ 1)
    b_boot[i,15] = est[5]
    b_boot[i,16] = (est[5]./10).^(est[12]+1)
    b_boot[i,17] = est[6]
    b_boot[i,18] = est[16]
    b_boot[i,19] = est[17]
    b_boot[i,20] = est[18]
    b_boot[i,21] = est[19]
    b_boot[i,22] = est[14]
    b_boot[i,23] = est[20]
    b_boot[i,24] = est[20] .* est[21]
    b_boot[i,25] = 1.0 - est[22]
end
betasigma_std_boot = [std(b_boot[:,j])             for j=1:25]
betasigma_boot_25  = [quantile(b_boot[:,j], 0.025) for j=1:25]
betasigma_boot_5   = [quantile(b_boot[:,j], 0.05)  for j=1:25]
betasigma_boot_10  = [quantile(b_boot[:,j], 0.1)   for j=1:25]
betasigma_boot_90  = [quantile(b_boot[:,j], 0.9)   for j=1:25]
betasigma_boot_95  = [quantile(b_boot[:,j], 0.95)  for j=1:25]
betasigma_boot_975 = [quantile(b_boot[:,j], 0.975) for j=1:25]

# Bootstrap Welfare Statistics
boot_welfare = CSV.read("data/bootstrap_welfare.csv", DataFrame, header=false)
boot_welfare = unique(boot_welfare)
boot_welfare = boot_welfare[:, (end-8):end]

welfare_std_boot = [std(boot_welfare[:,j])             for j=1:9]
welfare_boot_25  = [quantile(boot_welfare[:,j], 0.025) for j=1:9]
welfare_boot_5   = [quantile(boot_welfare[:,j], 0.05)  for j=1:9]
welfare_boot_10  = [quantile(boot_welfare[:,j], 0.1)   for j=1:9]
welfare_boot_90  = [quantile(boot_welfare[:,j], 0.9)   for j=1:9]
welfare_boot_95  = [quantile(boot_welfare[:,j], 0.95)  for j=1:9]
welfare_boot_975 = [quantile(boot_welfare[:,j], 0.975) for j=1:9]

# Column definitions are same as rows in Summary201609.xlsx.
CSV.write("data/bootstrap_estimates_$(vint).csv",         Tables.table(b_boot))
CSV.write("data/bootstrap_welfare_estimates_$(vint).csv", boot_welfare)
