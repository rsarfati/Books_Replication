# Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>
using CSV, DataFrames, Dates, Distributions, MAT, Statistics

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
mode      = 2   # Choose between available modes: 1 or 2 (see description above)
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
gamma0vec = vcat(quantile.(Gamma(0.5, 20), 0.005:0.01:0.895), 28:2:60, 64:4:100)
deltavec  = vcat(exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)), 3:2:20)
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

if mode == 2
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
for i = 1#:N_bs
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
    bdata12["N"] = length(bdata12["p"])
    bdata12["M"] = length(bdata12["first"])

    # bp
    bmktsizebp[i,:] = bp["mktsize"][bs_ind]
    bfirstbp[i,:]   = bp["first"][bs_ind]
    bcdindexbp[i,:] = bp["cdindex"][bs_ind]
    bendbp          = cumsum(bmktsizebp[i,:])
    bstartbp = vcat(1, bendbp[1:end-1] .+ 1)

    bbp = Dict{String,Any}()
    bbp["cdid"] = collect(1:236)  # bp.cdid[bs_ind]
    for x in [:p, :obsweight, :numlist, :condition, :localint, :popular, :conditiondif]
        bbp[String(x)] = bp[String(x)][bs_ind]
    end

    bbp["first"]  = bstartbp'
    bbp["cdindex"] = bendbp'
    bbp["mktsize"] = bmktsizebp[i,:]'
    bbp["N"] = length(bbp["p"])
    bbp["M"] = length(bbp["first"])

    # Optimization
    if mode == 1
        x0 = true_estimates[7:20]
        objectivefun = @(x) objective(x, x0, vars["distpara0"], gamma0vec, deltavec, bdata12, bdata09, bbp)
        x00 = x0
        x0[[3, 7]] = []

        x, fval = fminsearch(objectivefun, x0, optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5))
        estimates = [i, x[1:2], x00[3], x[3:5], x00[7], x[6:(length(x00)-2)]]
        CSV.write("bootstrap_estimates.csv", estimates)

        # The code below are just to show what bootstrap_betasigma.csv was. It appears not used.
        # xx = estimates
        # betasigma5new = [i, distpara0...
        #     [xx(1) xx(2)/(1+xx(1)) xx(3) ...
        #     xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
        #     xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)]
        # dlmwrite('bootstrap_betasigma.csv',betasigma5new,'delimiter',',','-append')

    elseif mode == 2
        xinitial = boot[i,:]
        llh, newdistpara, fother, fWF = full_model(xinitial, distpara0, gamma0vec, deltavec, bdata12, bdata09, bbp) # ignore input fWF to avoid welfare computation.

        # The code below just to show what bootstrap_distpara.csv was.
        # It will not be used.
        # result_w = [i,xinitial,newdistpara]
        # dlmwrite('bootstrap_distpara.csv',result_w,'delimiter',',','-append')

        wel09    = mean(fWF.AveWF09)
        wel12    = mean(fWF.AveWF12)
        weloff   = mean(fWF.WFbp)
        result_w = [i, xinitial, newdistpara, wel09, wel12, weloff]

        CSV.write("bootstrap_welfare.csv", result_w)

        # standard_errors.m will process bootstrap_distpara.csv and
        # bootstrap_welfare.csv to generate the numbers in Summary201609.xls.

        # The code below shows what bootstrap_results.csv was. It appears unused.

        # wel09 = mean(fWF.AveWF09)
        # wel12 = mean(fWF.AveWF12)
        # weloff = mean(fWF.WFbp)
        # output = [i,llh- 8223*log(20), newdistpara,estimates,clocktime]

        # dlmwrite('bootstrap_results.csv',output,'delimiter',',','-append')
    end
end

## Lines adapted from Masao's standard_errors.m, can be run after Mode 2.
# betasigma std -- bootstrap
boot = CSV.read("bootstrap_welfare.csv", DataFrame, header=true)
boot = unique(boot; dims=1)
boot = boot[:, 2:15]

distpara = CSV.read("bootstrap_welfare.csv", DataFrame, header=true)
distpara = unique(distpara; dims=1)
distpara = distpara[:, 16:21]

b_boot = zeros(size(boot,1),25)
for i = 1:size(boot,1)
    distpara0 = distpara[i,:]
    xx = boot[i,:]
    xtemp[i,:] = xx
    est = [distpara0...
        [xx[1] xx[2]/(1+xx[1]) xx[3] ...
         xx[4]*10*xx[7]/10/9.5^(-xx[6]-1) xx[5]*10*xx[7]/10/8^(-xx[6]-1) ...
         xx[6:11].*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx[12:13] xx[14]]
    bh[1] = est[1]
    bh[2] = est[2]
    bh[3] = est[3]
    bh[4] = est[7]
    bh[5] = est[8]
    bh[6] = est[13]
    bh[7] = est[15]
    bh[8] = est[4]
    bh[9] = est[12]
    bh[10] = est[9]
    bh[11] = est[10] .* est[9]
    bh[12] = (est[10] .* est[9] ./ 10) .^ (est[12] .+ 1)
    bh[13] = est[11] .* est[9]
    bh[14] = (est[11] .* est[9] ./ 10) .^ (est[12] .+ 1)
    bh[15] = est[5]
    bh[16] = (est[5]./10).^(est[12]+1)
    bh[17] = est[6]
    bh[18] = est[16]
    bh[19] = est[17]
    bh[20] = est[18]
    bh[21] = est[19]
    bh[22] = est[14]
    bh[23] = est[20]
    bh[24] = est[20] .* est[21]
    bh[25] = 1.0 - est[22]
    b_boot[i,:] = bh
end
betasigma_std_boot = zeros(25)

for j = 1:25
    betasigma_std_boot[j] = std(b_boot[:,j])
    betasigma_boot_25[j]  = quantile(b_boot[:,j], 0.025)
    betasigma_boot_5[j]   = quantile(b_boot[:,j], 0.05)
    betasigma_boot_10[j]  = quantile(b_boot[:,j], 0.1)
    betasigma_boot_90[j]  = quantile(b_boot[:,j], 0.9)
    betasigma_boot_95[j]  = quantile(b_boot[:,j], 0.95)
    betasigma_boot_975[j] = quantile(b_boot[:,j], 0.975)
end

#welfare std -- bootstrap
boot_welfare = CSV.read("bootstrap_welfare.csv", DataFrame, header=true)
boot_welfare = unique(boot; dims=1)
#boot_welfare     = boot_welfare(unique(boot_welfare(:, 1),'rows'),:)
#boot_welfare     = unique(boot_welfare, 'rows')
boot_welfare     = boot_welfare[:, (end-8):end]
welfare_std_boot = zeros(9)

for j = 1:9
    welfare_std_boot[j] = std(boot_welfare[:,j])
    welfare_boot_25[j]  = quantile(boot_welfare[:,j], 0.025)
    welfare_boot_5[j]   = quantile(boot_welfare[:,j], 0.05)
    welfare_boot_10[j]  = quantile(boot_welfare[:,j], 0.1)
    welfare_boot_90[j]  = quantile(boot_welfare[:,j], 0.9)
    welfare_boot_95[j]  = quantile(boot_welfare[:,j], 0.95)
    welfare_boot_975[j] = quantile(boot_welfare[:,j], 0.975)
end

# transformed bootstrap results are saved in two spreadsheets for analysis
# outside of matlab. Column definitions are same as rows in Summary201609.xlsx.
CSV.write("bootstrap_estimates_$(vint).csv",         b_boot)
CSV.write("bootstrap_welfare_estimates_$(vint).csv", boot_welfare)
