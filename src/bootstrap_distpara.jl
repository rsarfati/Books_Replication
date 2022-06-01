# Based on file <bootstrap_distpara_obtain_documented_Jan2018.m0>

## Input and output
# This file takes true data and randomly generated title-level index as
# inputs, contructs bootstrap dataset and runs the estimation, and output
# bootstrap_welfare.csv

# standard_errors.m (by Masao) will process this file to generate the
# bootstrap numbers in Summary201609.xls. I adapt the lines from that file
# at the end of this file so there is no need to look up that file.


## Estimates and welfare from true data

true_estimates = [0.32135	5.6459	14.855	1.1614	0.6486	1.9196	14.771
                  -2.4895	1	0.44004	0.32415	0.87235	0.5	0.25921	-9.1217
                  80.267	-13.647 1.7296	8.8188	0.92622	4.283	4.9097	0
                  7.8609	7.739	0.011111	6.8386	6.5027	0.028621]

# this is estimates and welfare from true data. It is copied from last row
# of results2.csv created by Masao. Except for the 5th and 6th elements,
# E[gamma_i] 2009 offline and betalocal, which according to my validation
# mode, is different from Masao's number.

## translating true_estimates to the "Estimates" column in Summary201609

#   # true_estimates contains three parts:
#     external_params = true_estimates(7:20) # parameters as first input to fullmodelllhWFAug22newtest2015_all, that is optimized with fullmodelllhWFAug22newtest2015_all being the objective function.
#     internal_params = true_estimates(7:20) # parameters as second input to fullmodelllhWFAug22newtest2015_all, that is optimized within fullmodelllhWFAug22newtest2015_all.
#     welfare_params = true_estimates(21:29) # welfare numbers computed by fullmodelllhWFAug22newtest2015_all.
#   # translate external_params and internal_params to row 2-26 in
#   # Summary201609. welfare_params' goes directly to row 28-36
#     xx = external_params
#     distpara0 = internal_params
#     est = [distpara0...
#         [xx(1) xx(2)/(1+xx(1)) xx(3) ...
#         xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
#         xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)]
#     bh(1) = est(1)
#     bh(2) = est(2)
#     bh(3) = est(3)
#     bh(4) = est(7)
#     bh(5) = est(8)
#     bh(6) = est(13)
#     bh(7) = est(15)
#     bh(8) = est(4)
#     bh(9) = est(12)
#     bh(10) = est(9)
#     bh(11) = est(10).*est(9)
#     bh(12) = (est(10).*est(9)./10).^(est(12)+1)
#     bh(13) = est(11).*est(9)
#     bh(14) = (est(11).*est(9)./10).^(est(12)+1)
#     bh(15) = est(5)
#     bh(16) = (est(5)./10).^(est(12)+1)
#     bh(17) = est(6)
#     bh(18) = est(16)
#     bh(19) = est(17)
#     bh(20) = est(18)
#     bh(21) = est(19)
#     bh(22) = est(14)
#     bh(23) = est(20)
#     bh(24) = est(20).*est(21)
#     bh(25) = 1-est(22)
#     # row 2-26 in Summary201609 is bh'


## Two modes
# this file has two modes. Mode 2 can be run after Mode 1.

# In the first mode, it cold starts to estimates parameters in the first
# input of the likelihood function fullmodelllhWFAug22newtest2015_all with
# bootstrap sample, and save the results in bootstrap_estimates.csv. These
# parameters are/transform to alpha, beta, gammaishape, gammaimean09,
# gammaimean12, eta, r, olp, c, lamda1, lamda2, betacond, betapop, olm,
# oltheta, naturaldisappear. This mode is very computationally intensive
# (involves optimization) and has to be run on servers and still takes
# weeks.

# In the second mode, it figures out the six parameters and welfare numbers
# that fullmodelllhWFAug22newtest2015_all can compute internally, reading
# the estimates from the first mode as input. This mode is quite fast
# without the welfare part, and can be run on a desktop in a couple hours.
# But it will take a day on servers if welfare output (WF) from
# fullmodelllhWFAug22newtest2015_all is included.

# As of Dec 2017, Mode 1 is done with the issue of duplicated book
# index in bootstrap sample, Mode 2 is done with the issue
# corrected, albeit taking Mode's 1's results as input.

current_mode = 2 # Alternatives: (1, 2)

## Parallel computing setup
# The following 3 lines are for reduced-form server.
# Can change resource parameters (walltime and mem) if necessary
# dcluster=parcluster
# dcluster.ResourceTemplate='-l nodes=^N^,software=MATLAB_Distrib_Comp_Engine+^N^,walltime=80:00:00,mem=64gb'
# dcluster.saveProfile

# This line is common for all servers or desktop. Choose the number accordingly.
#parpool(72)

## Load data
# Mat file created by Masao contains original data, random index for bootstrap.
load('DataToRun_pop09_boot.mat')

## Data Renaming and Setup pre-Bootstrap

gamma0vec = [gaminv(0.005:0.01:0.895,0.5,20) 28:2:60 64:4:100]
deltavec = [exp(norminv(0.01:0.02:0.91,-2,2)) 3:2:20]
data12 = data12nopop
data09 = data09nopop
bp = bpnopop

bmktsize09 = bootindex
bfirst09 = bootindex
bcdindex09 = bootindex

bmktsize12 = bootindex
bfirst12 = bootindex
bcdindex12 = bootindex

bmktsizebp = bootindex
bfirstbp = bootindex
bcdindexbp = bootindex

if current_mode == 2
    boot = csvread('bootstrap_estimates.csv')
    boot = unique(boot,'rows') # because mode 1 can be run on several servers simultaneously to save time, this line removes bootstrap runs that were duplicated on multiple servers.
    boot = boot(:,2:end)
end

## A validation mode:
# if one wants to reproduce the true data estimates or validate the
# estimation method is the same as we used for the true data estimates, the
# following tricks can be used to reproduce the true data set.

# in mode 1:
# bootindex(1,:) = 1:236
# in mode 2, in addition to the line above:
# boot(1,:) = results2(end,8:21)

## Bootstrap runs
beginning = 1
ending = 200

for i = beginning:ending
    ##  This block generates the bootstrap data.
    # 09
    bmktsize09[i,:] = data09.mktsize(bootindex[i,:])
    bfirst09[i,:] = data09.first(bootindex[i,:])
    bcdindex09[i,:] = data09.cdindex(bootindex[i,:])

    bdata09.p = zeros(sum(bmktsize09[i,:]),1)
    bdata09.cdid = zeros(sum(bmktsize09[i,:]),1)
    bdata09.obsweight = zeros(sum(bmktsize09[i,:]),1)
    bdata09.numlist = zeros(sum(bmktsize09[i,:]),1)
    bdata09.condition = zeros(sum(bmktsize09[i,:]),1)
    bdata09.localint = zeros(sum(bmktsize09[i,:]),1)
    bdata09.popular = zeros(sum(bmktsize09[i,:]),1)
    bdata09.pdif = zeros(sum(bmktsize09[i,:]),1)
    bdata09.conditiondif = zeros(sum(bmktsize09[i,:]),1)
    bdata09.basecond = zeros(sum(bmktsize09[i,:]),1)


    bend09 = cumsum(bmktsize09[i,:])
    bstart09temp = bend09+1
    bstart09 = [1 bstart09temp(1:(length(bstart09temp)-1))]
    for j = 1:length(bfirst09[i,:])
        bdata09.p(bstart09[j]:bend09[j]) = data09.p(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.cdid(bstart09[j]:bend09[j]) = repmat(j, bend09[j] - bstart09[j] + 1,1)
        bdata09.obsweight(bstart09[j]:bend09[j]) = data09.obsweight(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.numlist(bstart09[j]:bend09[j]) = data09.numlist(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.condition(bstart09[j]:bend09[j]) = data09.condition(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.localint(bstart09[j]:bend09[j]) = data09.localint(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.popular(bstart09[j]:bend09[j]) = data09.popular(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.pdif(bstart09[j]:bend09[j]) = data09.pdif(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.conditiondif(bstart09[j]:bend09[j]) = data09.conditiondif(bfirst09[i,j]:bcdindex09[i,j])
        bdata09.basecond(bstart09[j]:bend09[j]) = data09.basecond(bfirst09[i,j]:bcdindex09[i,j])
    end
    bdata09.first = bstart09'
    bdata09.cdindex = bend09'
    bdata09.mktsize = bmktsize09[i,:]'
    bdata09.N = length(bdata09.p)
    bdata09.M = length(bdata09.first)


    # 12
    bmktsize12[i,:] = data12.mktsize(bootindex[i,:])
    bfirst12[i,:] = data12.first(bootindex[i,:])
    bcdindex12[i,:] = data12.cdindex(bootindex[i,:])

    bdata12.p = zeros(sum(bmktsize12[i,:]),1)
    bdata12.cdid = zeros(sum(bmktsize12[i,:]),1)
    bdata12.obsweight = zeros(sum(bmktsize12[i,:]),1)
    bdata12.numlist = zeros(sum(bmktsize12[i,:]),1)
    bdata12.condition = zeros(sum(bmktsize12[i,:]),1)
    bdata12.localint = zeros(sum(bmktsize12[i,:]),1)
    bdata12.popular = zeros(sum(bmktsize12[i,:]),1)
    bdata12.pdif = zeros(sum(bmktsize12[i,:]),1)
    bdata12.conditiondif = zeros(sum(bmktsize12[i,:]),1)
    bdata12.basecond = zeros(sum(bmktsize12[i,:]),1)
    bdata12.disappear = zeros(sum(bmktsize12[i,:]),1)

    bend12 = cumsum(bmktsize12[i,:])
    bstart12temp = bend12+1
    bstart12 = [1 bstart12temp(1:(length(bstart12temp)-1))]
    for j = 1:length(bfirst12[i,:])
        bdata12.p(bstart12[j]:bend12[j]) = data12.p(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.cdid(bstart12[j]:bend12[j]) = repmat(j, bend12[j] - bstart12[j] + 1,1)
        bdata12.obsweight(bstart12[j]:bend12[j]) = data12.obsweight(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.numlist(bstart12[j]:bend12[j]) = data12.numlist(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.condition(bstart12[j]:bend12[j]) = data12.condition(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.localint(bstart12[j]:bend12[j]) = data12.localint(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.popular(bstart12[j]:bend12[j]) = data12.popular(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.pdif(bstart12[j]:bend12[j]) = data12.pdif(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.conditiondif(bstart12[j]:bend12[j]) = data12.conditiondif(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.basecond(bstart12[j]:bend12[j]) = data12.basecond(bfirst12[i,j]:bcdindex12[i,j])
        bdata12.disappear(bstart12[j]:bend12[j]) = data12.disappear(bfirst12[i,j]:bcdindex12[i,j])
    end
    bdata12.first = bstart12'
    bdata12.cdindex = bend12'
    bdata12.mktsize = bmktsize12[i,:]'
    bdata12.N = length(bdata12.p)
    bdata12.M = length(bdata12.first)

    # bp
    bmktsizebp[i,:] = bp.mktsize(bootindex[i,:])
    bfirstbp[i,:] = bp.first(bootindex[i,:])
    bcdindexbp[i,:] = bp.cdindex(bootindex[i,:])
    bendbp = cumsum(bmktsizebp[i,:])
    bstartbptemp = bendbp+1
    bstartbp = [1 bstartbptemp(1:(length(bstartbptemp)-1))]

    bbp.p = bp.p(bootindex[i,:])
    bbp.cdid = (1:236)'  # bp.cdid(bootindex[i,:])
    bbp.obsweight = bp.obsweight(bootindex[i,:])
    bbp.numlist = bp.numlist(bootindex[i,:])
    bbp.condition = bp.condition(bootindex[i,:])
    bbp.localint = bp.localint(bootindex[i,:])
    bbp.popular = bp.popular(bootindex[i,:])
    bbp.conditiondif = bp.conditiondif(bootindex[i,:])

    bbp.first = bstartbp'
    bbp.cdindex = bendbp'
    bbp.mktsize = bmktsizebp[i,:]'
    bbp.N = length(bbp.p)
    bbp.M = length(bbp.first)

    ## This part execute optimization/computation
    if current_mode == 1
        x0 = true_estimates[7:20]
        objectivefun = @(x) objective(x,x0,distpara0,gamma0vec,deltavec,bdata12,bdata09,bbp)
        x00 = x0
        x0[[3, 7]] = []
        [x,fval] = fminsearch(objectivefun,x0,optimset('MaxFunEvals',1e5,'MaxIter',1e5))
        estimates = [i,x[1:2] x00[3] x[3:5] x00[7] x[6:(length(x00)-2)]]
        dlmwrite('bootstrap_estimates.csv',estimates,'delimiter',',','-append')

        # The code below are just to show what bootstrap_betasigma.csv was. It appears not used.
        # xx = estimates
        # betasigma5new = [i, distpara0...
        #     [xx(1) xx(2)/(1+xx(1)) xx(3) ...
        #     xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
        #     xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)]
        # dlmwrite('bootstrap_betasigma.csv',betasigma5new,'delimiter',',','-append')

    elseif current_mode == 2
        xinitial = boot[i,:]
        [llh, newdistpara,fother,fWF] = fullmodelllhWFAug22newtest2015(xinitial,distpara0,gamma0vec,deltavec, bdata12,bdata09,bbp) # ignore input fWF to avoid welfare computation.

        #     The code below are just to show what bootstrap_distpara.csv was.
        #     It will not be used.
        # result_w = [i,xinitial,newdistpara]
        #
        # dlmwrite('bootstrap_distpara.csv',result_w,'delimiter',',','-append')

        wel09 = mean(fWF.AveWF09)
        wel12 = mean(fWF.AveWF12)
        weloff = mean(fWF.WFbp)
        result_w = [i,xinitial,newdistpara,wel09,wel12,weloff]


        dlmwrite('bootstrap_welfare.csv',result_w,'delimiter',',','-append')

        # standard_errors.m will process bootstrap_distpara.csv and
        # bootstrap_welfare.csv to generate the numbers in Summary201609.xls.
        #
        # The code below are just to show what bootstrap_results.csv was. It appears not used.
        # clocktime = clock
        # clocktime = clocktime(1:5)
        #
        # wel09 = mean(fWF.AveWF09)
        # wel12 = mean(fWF.AveWF12)
        # weloff = mean(fWF.WFbp)
        # output = [i,llh- 8223*log(20), newdistpara,estimates,clocktime]
        #
        # dlmwrite('bootstrap_results.csv',output,'delimiter',',','-append')
    end
end

## Lines adapted from Masao's standard_errors.m, can be run after Mode 2.
# betasigma std -- bootstrap
boot = csvread('bootstrap_welfare.csv')
boot = unique(boot,'rows')
boot = boot(:,2:15)
distpara = csvread('bootstrap_welfare.csv')
distpara = unique(distpara,'rows')
distpara = distpara(:,16:21)

b_boot = zeros(size(boot,1),25)
for i = 1:size(boot,1)
    distpara0 = distpara[i,:]
    xx = boot[i,:]
    xtemp[i,:] = xx
    est = [distpara0...
        [xx(1) xx(2)/(1+xx(1)) xx(3) ...
        xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
        xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)]
    bh(1) = est(1)
    bh(2) = est(2)
    bh(3) = est(3)
    bh(4) = est(7)
    bh(5) = est(8)
    bh(6) = est(13)
    bh(7) = est(15)
    bh(8) = est(4)
    bh(9) = est(12)
    bh(10) = est(9)
    bh(11) = est(10).*est(9)
    bh(12) = (est(10).*est(9)./10).^(est(12)+1)
    bh(13) = est(11).*est(9)
    bh(14) = (est(11).*est(9)./10).^(est(12)+1)
    bh(15) = est(5)
    bh(16) = (est(5)./10).^(est(12)+1)
    bh(17) = est(6)
    bh(18) = est(16)
    bh(19) = est(17)
    bh(20) = est(18)
    bh(21) = est(19)
    bh(22) = est(14)
    bh(23) = est(20)
    bh(24) = est(20).*est(21)
    bh(25) = 1-est(22)
    b_boot[i,:] = bh
end
betasigma_std_boot = zeros(25,1)
for j = 1:25
    betasigma_std_boot[j] = std(b_boot(:,j))
    betasigma_boot_25[j] = quantile(b_boot(:,j),0.025)
    betasigma_boot_5[j] = quantile(b_boot(:,j),0.05)
    betasigma_boot_10[j] = quantile(b_boot(:,j),0.1)
    betasigma_boot_90[j] = quantile(b_boot(:,j),0.9)
    betasigma_boot_95[j] = quantile(b_boot(:,j),0.95)
    betasigma_boot_975[j] = quantile(b_boot(:,j),0.975)
end

#welfare std -- bootstrap
boot_welfare = csvread('bootstrap_welfare.csv')
boot_welfare = boot_welfare(unique(boot_welfare(:,1),'rows'),:)
boot_welfare = unique(boot_welfare,'rows')
boot_welfare = boot_welfare(:,(end-8):end)
welfare_std_boot = zeros(9,1)
for j = 1:9
    welfare_std_boot[j] = std(boot_welfare(:,j))
    welfare_boot_25[j] = quantile(boot_welfare(:,j),0.025)
    welfare_boot_5[j] = quantile(boot_welfare(:,j),0.05)
    welfare_boot_10[j] = quantile(boot_welfare(:,j),0.1)
    welfare_boot_90[j] = quantile(boot_welfare(:,j),0.9)
    welfare_boot_95[j] = quantile(boot_welfare(:,j),0.95)
    welfare_boot_975[j] = quantile(boot_welfare(:,j),0.975)
end

# transformed bootstrap results are saved in two spreadsheets for analysis
# outside of matlab. Column definitions are same as rows in Summary201609.xlsx.
csvwrite('bootstrap_estimates2017_1105.csv',b_boot)
csvwrite('bootstrap_welfare_estimates2017_1105.csv',boot_welfare)
