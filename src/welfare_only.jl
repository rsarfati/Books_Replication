## Welfare only
load('DataToRun_pop09.mat')
betasigma5 = [ distpara0...
    [x0(1) x0(2)/(1+x0(1)) x0(3) ...
    x0(4)*10*x0(7)/10/9.5^(-x0(6)-1) x0(5)*10*x0(7)/10/8^(-x0(6)-1) ...
    x0(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] x0(12:13) x0(14)];
gamma0vec = [gaminv(0.005:0.01:0.895,0.5,20) 28:2:60 64:4:100];
deltavec = [exp(norminv(0.01:0.02:0.91,-2,2)) 3:2:20];
data12 = data12nopop;
data09 = data09nopop;
bp = bpnopop;


%%
parpool(20)

estimates = csvread('estimates.csv');
estimates = estimates(1,:);

[llh, newdistpara,fother,fWF] = fullmodelllhWFAug22newtest2015(estimates,distpara0,gamma0vec,deltavec, data12,data09,bp);
clocktime = clock;
clocktime = clocktime(1:5);

save('fWF.mat','fWF')

wel09 = mean(fWF.AveWF09);
wel12 = mean(fWF.AveWF12);
weloff = mean(fWF.WFbp);
output = [llh- 8223*log(20), newdistpara,estimates,weloff,wel09,wel12,clocktime];

dlmwrite('results2.csv',output,'delimiter',',','-append');


% x0 and newdistpara can be transformed to the estimated paramters in the
% 10/15/2015 column of resultssummary.xlsx. (some elements are already the
% parameters reported in the spreadsheet, others needs some simple
% transformation)


% llh - 8223*log(20) is the likelihood in the 10/15/2015 column of resultssummary.xlsx
