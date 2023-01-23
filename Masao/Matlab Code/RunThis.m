clear
clear global
close all
load('DataToRun_pop09.mat')
est = csvread('Bootstrapp_data/bootstrap_estimates.csv');
%x0 = est(end,:);
%x0 = mean(est(:,2:end),1);
x0  = [4.519934831
    -2.260332606
    1
    0.491626599
    0.19667731
    0.866176021
    0.5
    0.275033092
    -12.74106448
    95.40782168
    -19.20384749
    1.68921498
    9.458397134
    0.931646012]';

gamma0vec = [gaminv(0.005:0.01:0.895,0.5,20) 28:2:60 64:4:100];
deltavec = [exp(norminv(0.01:0.02:0.91,-2,2)) 3:2:20];
data12 = data12nopop;
data09 = data09nopop;
bp = bpnopop;


%%
%parpool(20)
%%

% objectivefun = @(x) objective(x,x0,distpara0,gamma0vec,deltavec,data12,data09,bp);
% x00 = x0;
% 
% x0([3 7]) = [];
% 
% [x,fval] = fminsearch(objectivefun,x0,optimset('Display','iter','MaxFunEvals',1e4,'MaxIter',1e4));


%%
%estimates = [x(1:2) x00(3) x(3:5) x00(7) x(6:(length(x00)-2))];
%dlmwrite('estimates.csv',estimates,'delimiter',',','-append');
%xx = estimates;
%betasigma5new = [ distpara0...
%     [xx(1) xx(2)/(1+xx(1)) xx(3) ...
%     xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
%     xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
% dlmwrite('betasigma_record.csv',betasigma5new,'delimiter',',','-append');
% 
[llh, newdistpara, fother, fWF] = fullmodelllhWFAug22newtest2015(x0,distpara0,gamma0vec,deltavec,data12,data09,bp);
clocktime = clock;
clocktime = clocktime(1:5);

wel09 = mean(fWF.AveWF09);
wel12 = mean(fWF.AveWF12);
weloff = mean(fWF.WFbp);
output = [llh- 8223*log(20), newdistpara, estimates, weloff, wel09, wel12, clocktime];

dlmwrite('results2.csv',output,'delimiter',',','-append');
%save('fWF.mat','fWF')
% x0 and newdistpara can be transformed to the estimated paramters in the
% 10/15/2015 column of resultssummary.xlsx. (some elements are already the
% parameters reported in the spreadsheet, others needs some simple
% transformation)


% llh - 8223*log(20) is the likelihood in the 10/15/2015 column of resultssummary.xlsx
