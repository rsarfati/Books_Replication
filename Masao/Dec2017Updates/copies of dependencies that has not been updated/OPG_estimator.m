clear
clear global
close all
%parpool(20)
load('DataToRun_pop09.mat')
est = csvread('estimates.csv');
x0 = est(end,:);
distpara = csvread('estimates_distpara.csv');
distpara0 = distpara(end,:);
gamma0vec = [gaminv(0.005:0.01:0.895,0.5,20) 28:2:60 64:4:100];
deltavec = [exp(norminv(0.01:0.02:0.91,-2,2)) 3:2:20];
data12 = data12nopop;
data09 = data09nopop;
bp = bpnopop;


M = bp.M;
llh_single = zeros(M,1);
step = 0.0001;

parameter = [x0, distpara0];
parameter(3) = [];
parameter(6) = [];
npara = length(parameter);
nx = length(x0)-2;
nd = length(distpara0);
%% start loop
llh_side = zeros(2,1);
dllh = zeros(M,npara);
V = zeros(npara,npara);
for k = 1:M
    for np = 1:npara
        for j = 1:2
            dparam = parameter;
            if j == 1
                dparam(np) = dparam(np) + step;
            else
                dparam(np) = dparam(np) - step;
            end
            llh_side(j) = fullmodelllhWFAug22newtest2015_single([dparam(1:2),x0(3),dparam(3:5),x0(7),dparam(6:nx)],dparam((nx+1):end),gamma0vec,deltavec, data12,data09,bp,k);
        end
        dllh(k,np) = (llh_side(1)-llh_side(2))./(2*step);
    end
    V = V + dllh(k,:)'*dllh(k,:);
end
OPGmatrix = inv(1./M.*V);
std = sqrt(diag(OPGmatrix));
%% write standard error
dlmwrite('OPGestimator.csv',std,'delimiter',',','-append');


