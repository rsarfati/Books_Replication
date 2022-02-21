%% Initial Draw
% clear
% load('DataToRun_pop09.mat')
% data12 = data12nopop;
% data09 = data09nopop;
% bp = bpnopop;
%
% id = unique(data09.cdid);
% rng(10);
% Nboot = 200;
% bootindex = zeros(Nboot,length(id));
% bcdid = bootindex;
%
% for i = 1:Nboot
%     bootindex(i,:) = datasample(1:length(id),length(id));
%     bcdid(i,:) = id(bootindex(i,:));
% end
% clear data12 data09 id bp
% save('DataToRun_pop09_boot.mat')
%% Initial setting
clear
clear global
close all

parpool(20);
jj = 6;
NN = 30;


load('DataToRun_pop09_boot.mat')
est = csvread('estimates.csv');
xinitial = est(end,:);
gamma0vec = [gaminv(0.005:0.01:0.895,0.5,20) 28:2:60 64:4:100];
deltavec = [exp(norminv(0.01:0.02:0.91,-2,2)) 3:2:20];
data12 = data12nopop;
data09 = data09nopop;
bp = bpnopop;



%%


bmktsize09 = bootindex;
bfirst09 = bootindex;
bcdindex09 = bootindex;

bmktsize12 = bootindex;
bfirst12 = bootindex;
bcdindex12 = bootindex;

bmktsizebp = bootindex;
bfirstbp = bootindex;
bcdindexbp = bootindex;

%beginning = ((jj-1)*NN+1);
%ending = (jj*NN);
beginning = 140;
ending = 140;
for i = beginning:ending
    done = csvread('bootstrap_estimates.csv');
    donelist = done(:,1);
    if any(abs(i-donelist)<0.01) == 0
        
        % 09
        bmktsize09(i,:) = data09.mktsize(bootindex(i,:));
        bfirst09(i,:) = data09.first(bootindex(i,:));
        bcdindex09(i,:) = data09.cdindex(bootindex(i,:));
        
        bdata09.p = zeros(sum(bmktsize09(i,:)),1);
        bdata09.cdid = zeros(sum(bmktsize09(i,:)),1);
        bdata09.obsweight = zeros(sum(bmktsize09(i,:)),1);
        bdata09.numlist = zeros(sum(bmktsize09(i,:)),1);
        bdata09.condition = zeros(sum(bmktsize09(i,:)),1);
        bdata09.localint = zeros(sum(bmktsize09(i,:)),1);
        bdata09.popular = zeros(sum(bmktsize09(i,:)),1);
        bdata09.pdif = zeros(sum(bmktsize09(i,:)),1);
        bdata09.conditiondif = zeros(sum(bmktsize09(i,:)),1);
        bdata09.basecond = zeros(sum(bmktsize09(i,:)),1);
        
        
        bend09 = cumsum(bmktsize09(i,:));
        bstart09temp = bend09+1;
        bstart09 = [1 bstart09temp(1:(length(bstart09temp)-1))];
        for j = 1:length(bfirst09(i,:))
            bdata09.p(bstart09(j):bend09(j)) = data09.p(bfirst09(i,j):bcdindex09(i,j));
            bdata09.cdid(bstart09(j):bend09(j)) = data09.cdid(bfirst09(i,j):bcdindex09(i,j));
            bdata09.obsweight(bstart09(j):bend09(j)) = data09.obsweight(bfirst09(i,j):bcdindex09(i,j));
            bdata09.numlist(bstart09(j):bend09(j)) = data09.numlist(bfirst09(i,j):bcdindex09(i,j));
            bdata09.condition(bstart09(j):bend09(j)) = data09.condition(bfirst09(i,j):bcdindex09(i,j));
            bdata09.localint(bstart09(j):bend09(j)) = data09.localint(bfirst09(i,j):bcdindex09(i,j));
            bdata09.popular(bstart09(j):bend09(j)) = data09.popular(bfirst09(i,j):bcdindex09(i,j));
            bdata09.pdif(bstart09(j):bend09(j)) = data09.pdif(bfirst09(i,j):bcdindex09(i,j));
            bdata09.conditiondif(bstart09(j):bend09(j)) = data09.conditiondif(bfirst09(i,j):bcdindex09(i,j));
            bdata09.basecond(bstart09(j):bend09(j)) = data09.basecond(bfirst09(i,j):bcdindex09(i,j));
        end
        bdata09.first = bstart09';
        bdata09.cdindex = bend09';
        bdata09.mktsize = bmktsize09(i,:)';
        bdata09.N = length(bdata09.p);
        bdata09.M = length(bdata09.first);
        
        
        % 12
        bmktsize12(i,:) = data12.mktsize(bootindex(i,:));
        bfirst12(i,:) = data12.first(bootindex(i,:));
        bcdindex12(i,:) = data12.cdindex(bootindex(i,:));
        
        bdata12.p = zeros(sum(bmktsize12(i,:)),1);
        bdata12.cdid = zeros(sum(bmktsize12(i,:)),1);
        bdata12.obsweight = zeros(sum(bmktsize12(i,:)),1);
        bdata12.numlist = zeros(sum(bmktsize12(i,:)),1);
        bdata12.condition = zeros(sum(bmktsize12(i,:)),1);
        bdata12.localint = zeros(sum(bmktsize12(i,:)),1);
        bdata12.popular = zeros(sum(bmktsize12(i,:)),1);
        bdata12.pdif = zeros(sum(bmktsize12(i,:)),1);
        bdata12.conditiondif = zeros(sum(bmktsize12(i,:)),1);
        bdata12.basecond = zeros(sum(bmktsize12(i,:)),1);
        bdata12.disappear = zeros(sum(bmktsize12(i,:)),1);
        
        bend12 = cumsum(bmktsize12(i,:));
        bstart12temp = bend12+1;
        bstart12 = [1 bstart12temp(1:(length(bstart12temp)-1))];
        for j = 1:length(bfirst12(i,:))
            bdata12.p(bstart12(j):bend12(j)) = data12.p(bfirst12(i,j):bcdindex12(i,j));
            bdata12.cdid(bstart12(j):bend12(j)) = data12.cdid(bfirst12(i,j):bcdindex12(i,j));
            bdata12.obsweight(bstart12(j):bend12(j)) = data12.obsweight(bfirst12(i,j):bcdindex12(i,j));
            bdata12.numlist(bstart12(j):bend12(j)) = data12.numlist(bfirst12(i,j):bcdindex12(i,j));
            bdata12.condition(bstart12(j):bend12(j)) = data12.condition(bfirst12(i,j):bcdindex12(i,j));
            bdata12.localint(bstart12(j):bend12(j)) = data12.localint(bfirst12(i,j):bcdindex12(i,j));
            bdata12.popular(bstart12(j):bend12(j)) = data12.popular(bfirst12(i,j):bcdindex12(i,j));
            bdata12.pdif(bstart12(j):bend12(j)) = data12.pdif(bfirst12(i,j):bcdindex12(i,j));
            bdata12.conditiondif(bstart12(j):bend12(j)) = data12.conditiondif(bfirst12(i,j):bcdindex12(i,j));
            bdata12.basecond(bstart12(j):bend12(j)) = data12.basecond(bfirst12(i,j):bcdindex12(i,j));
            bdata12.disappear(bstart12(j):bend12(j)) = data12.disappear(bfirst12(i,j):bcdindex12(i,j));
        end
        bdata12.first = bstart12';
        bdata12.cdindex = bend12';
        bdata12.mktsize = bmktsize12(i,:)';
        bdata12.N = length(bdata12.p);
        bdata12.M = length(bdata12.first);
        
        % bp
        bmktsizebp(i,:) = bp.mktsize(bootindex(i,:));
        bfirstbp(i,:) = bp.first(bootindex(i,:));
        bcdindexbp(i,:) = bp.cdindex(bootindex(i,:));
        bendbp = cumsum(bmktsizebp(i,:));
        bstartbptemp = bendbp+1;
        bstartbp = [1 bstartbptemp(1:(length(bstartbptemp)-1))];
        
        bbp.p = bp.p(bootindex(i,:));
        bbp.cdid = bp.cdid(bootindex(i,:));
        bbp.obsweight = bp.obsweight(bootindex(i,:));
        bbp.numlist = bp.numlist(bootindex(i,:));
        bbp.condition = bp.condition(bootindex(i,:));
        bbp.localint = bp.localint(bootindex(i,:));
        bbp.popular = bp.popular(bootindex(i,:));
        bbp.conditiondif = bp.conditiondif(bootindex(i,:));
        
        bbp.first = bstartbp';
        bbp.cdindex = bendbp';
        bbp.mktsize = bmktsizebp(i,:)';
        bbp.N = length(bbp.p);
        bbp.M = length(bbp.first);
        
        
        % estimatation
        x0 = xinitial;
        objectivefun = @(x) objective(x,x0,distpara0,gamma0vec,deltavec, bdata12,bdata09,bbp);
        x00 = x0;
        
        x0([3 7]) = [];
        
        [x,fval] = fminsearch(objectivefun,x0,optimset('MaxFunEvals',1e5,'MaxIter',1e5));
        estimates = [i,x(1:2) x00(3) x(3:5) x00(7) x(6:(length(x00)-2))];
        dlmwrite('bootstrap_estimates.csv',estimates,'delimiter',',','-append');
        xx = estimates;
        betasigma5new = [i, distpara0...
            [xx(1) xx(2)/(1+xx(1)) xx(3) ...
            xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
            xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
        dlmwrite('bootstrap_betasigma_.csv',betasigma5new,'delimiter',',','-append');
        
        %[llh, newdistpara,fother,fWF] = fullmodelllhWFAug22newtest2015(estimates,distpara0,gamma0vec,deltavec, bdata12,bdata09,bbp);
        %clocktime = clock;
        %clocktime = clocktime(1:5);
        
        %wel09 = mean(fWF.AveWF09);
        %wel12 = mean(fWF.AveWF12);
        %weloff = mean(fWF.WFbp);
        %output = [i,llh- 8223*log(20), newdistpara,estimates,clocktime];
        
        %dlmwrite('bootstrap_results.csv',output,'delimiter',',','-append');
        
        
    end
end
