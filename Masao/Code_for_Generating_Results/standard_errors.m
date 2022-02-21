clear
close all
clear global
load('DataToRun_pop09.mat')

boot = csvread('bootstrap_distpara.csv');
boot = boot(:,2:end);
np = size(boot,2);
std_boot = zeros(np,1);
for i = 1:np
    std_boot(i) = std(boot(:,i));
end


M = bpnopop.M;
V = csvread('OPG_matrix.csv');
V = 1/M*V;
std_opgtemp = sqrt(diag(V));

std_opg = zeros(np,1);
std_opg(1:2) = std_opgtemp(1:2);
std_opg(3) = 0;
std_opg(4:6) = std_opgtemp(3:5);
std_opg(8:end) = std_opgtemp(6:end);

std_list = [std_opg,std_boot];


estimates = csvread('results2.csv');
x(1:14) = estimates(end,8:21);
x(15:20) = estimates(end,2:7);

result = [x',std_list];

% dlmwrite([pwd '/result/raw_estimates.csv'],result,'delimiter',',');


%% in terms of betasigma --- estimates
distpara0 = x(15:20);
xx = x(1:14);
est = [distpara0...
    [xx(1) xx(2)/(1+xx(1)) xx(3) ...
    xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
    xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
bh(1) = est(1);
bh(2) = est(2);
bh(3) = est(3);
bh(4) = est(7);
bh(5) = est(8);
bh(6) = est(13);
bh(7) = est(15);
bh(8) = est(4);
bh(9) = est(12);
bh(10) = est(9);
bh(11) = est(10).*est(9);
bh(12) = (est(10).*est(9)./10).^(est(12)+1);
bh(13) = est(11).*est(9);
bh(14) = (est(11).*est(9)./10).^(est(12)+1);
bh(15) = est(5);
bh(16) = (est(5)./10).^(est(12)+1);
bh(17) = est(6);
bh(18) = est(16);
bh(19) = est(17);
bh(20) = est(18);
bh(21) = est(19);
bh(22) = est(14);
bh(23) = est(20);
bh(24) = est(20).*est(21);
bh(25) = 1-est(22);

beta_sigma_estimate = bh;

%% betasigma std OPG
npara = length(x);
bh_temp = zeros(25,2);
d_bh = zeros(25,npara);
step = 0.0001;

for np = 1:npara
    for j = 1:2
        dparam = x;
        if j == 1
            dparam(np) = dparam(np) + step;
        else
            dparam(np) = dparam(np) - step;
        end
        
        distpara0 = dparam(15:20);
        xx = dparam(1:14);
        est = [distpara0...
            [xx(1) xx(2)/(1+xx(1)) xx(3) ...
            xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
            xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
        bh(1) = est(1);
        bh(2) = est(2);
        bh(3) = est(3);
        bh(4) = est(7);
        bh(5) = est(8);
        bh(6) = est(13);
        bh(7) = est(15);
        bh(8) = est(4);
        bh(9) = est(12);
        bh(10) = est(9);
        bh(11) = est(10).*est(9);
        bh(12) = (est(10).*est(9)./10).^(est(12)+1);
        bh(13) = est(11).*est(9);
        bh(14) = (est(11).*est(9)./10).^(est(12)+1);
        bh(15) = est(5);
        bh(16) = (est(5)./10).^(est(12)+1);
        bh(17) = est(6);
        bh(18) = est(16);
        bh(19) = est(17);
        bh(20) = est(18);
        bh(21) = est(19);
        bh(22) = est(14);
        bh(23) = est(20);
        bh(24) = est(20).*est(21);
        bh(25) = 1-est(22);
        bh_temp(:,j) = bh;
    end
    d_bh(:,np) = (bh_temp(:,2) - bh_temp(:,1))./(2*step);
end
d_bh_new = d_bh;
d_bh_new(:,3) = [];
d_bh_new(:,6) = [];



betasigma_std_opg = zeros(25,1);
for j = 1:25
    betasigma_std_opg(j) = sqrt(d_bh_new(j,:)*V*d_bh_new(j,:)');
end


%% std for welfare - OPG
d_welfare = csvread('welfare_d.csv');
d_welfare = d_welfare';
welfare_std_opg = zeros(9,1);
for j = 1:9
    welfare_std_opg(j) = sqrt(d_welfare(j,:)*V*d_welfare(j,:)');
end


%% betasigma std -- bootstrap


boot = csvread('bootstrap_distpara.csv');
boot = unique(boot,'rows');
boot = boot(:,2:15);
distpara = csvread('bootstrap_distpara.csv');
distpara = unique(distpara,'rows');
distpara = distpara(:,16:end);

b_boot = zeros(size(boot,1),25);
for i = 1:size(boot,1)
    distpara0 = distpara(i,:);
    xx = boot(i,:);
    xtemp(i,:) = xx;
    est = [distpara0...
        [xx(1) xx(2)/(1+xx(1)) xx(3) ...
        xx(4)*10*xx(7)/10/9.5^(-xx(6)-1) xx(5)*10*xx(7)/10/8^(-xx(6)-1) ...
        xx(6:11).*[1 0.1 1 0.1 0.01 0.1 ] 0 0] xx(12:13) xx(14)];
    bh(1) = est(1);
    bh(2) = est(2);
    bh(3) = est(3);
    bh(4) = est(7);
    bh(5) = est(8);
    bh(6) = est(13);
    bh(7) = est(15);
    bh(8) = est(4);
    bh(9) = est(12);
    bh(10) = est(9);
    bh(11) = est(10).*est(9);
    bh(12) = (est(10).*est(9)./10).^(est(12)+1);
    bh(13) = est(11).*est(9);
    bh(14) = (est(11).*est(9)./10).^(est(12)+1);
    bh(15) = est(5);
    bh(16) = (est(5)./10).^(est(12)+1);
    bh(17) = est(6);
    bh(18) = est(16);
    bh(19) = est(17);
    bh(20) = est(18);
    bh(21) = est(19);
    bh(22) = est(14);
    bh(23) = est(20);
    bh(24) = est(20).*est(21);
    bh(25) = 1-est(22);
    b_boot(i,:) = bh;
end
betasigma_std_boot = zeros(25,1);
for j = 1:25
    betasigma_std_boot(j) = std(b_boot(:,j));
    betasigma_boot_25(j) = quantile(b_boot(:,j),0.025);
    betasigma_boot_5(j) = quantile(b_boot(:,j),0.05);
    betasigma_boot_10(j) = quantile(b_boot(:,j),0.1);
    betasigma_boot_90(j) = quantile(b_boot(:,j),0.9);
    betasigma_boot_95(j) = quantile(b_boot(:,j),0.95);
    betasigma_boot_975(j) = quantile(b_boot(:,j),0.975);
end

%% welfare std -- bootstrap
boot_welfare = csvread('bootstrap_welfare.csv');
boot_welfare = boot_welfare(unique(boot_welfare(:,1),'rows'),:);
boot_welfare = unique(boot_welfare,'rows');
boot_welfare = boot_welfare(:,(end-8):end);
welfare_std_boot = zeros(9,1);
for j = 1:9
    welfare_std_boot(j) = std(boot_welfare(:,j));
    welfare_boot_25(j) = quantile(boot_welfare(:,j),0.025);
    welfare_boot_5(j) = quantile(boot_welfare(:,j),0.05);
    welfare_boot_10(j) = quantile(boot_welfare(:,j),0.1);
    welfare_boot_90(j) = quantile(boot_welfare(:,j),0.9);
    welfare_boot_95(j) = quantile(boot_welfare(:,j),0.95);
    welfare_boot_975(j) = quantile(boot_welfare(:,j),0.975);
end



%% write
betasigma_summary = [beta_sigma_estimate' betasigma_std_opg betasigma_std_boot ...
    betasigma_boot_25'  betasigma_boot_975' betasigma_boot_5'  betasigma_boot_95' betasigma_boot_10'  betasigma_boot_90'];
% dlmwrite([pwd '/result/betasigma_estimates.csv'],betasigma_summary,'delimiter',',');
welfare_std= [welfare_std_opg welfare_std_boot welfare_boot_25' welfare_boot_975' welfare_boot_5' welfare_boot_95'...
     welfare_boot_10' welfare_boot_90'];
% dlmwrite([pwd '/result/welfare_std.csv'],welfare_std,'delimiter',',');
dlmwrite([pwd '/result/bootstrap_estimates2017_1105.csv'],b_boot,'delimiter',',');
dlmwrite([pwd '/result/bootstrap_welfare_estimates2017_1105.csv'],boot_welfare,'delimiter',',');