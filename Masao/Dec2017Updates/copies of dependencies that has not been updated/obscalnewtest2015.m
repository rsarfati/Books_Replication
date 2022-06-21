function [lip,gamma2,gamma1,gamma0,D0,Dm,pi,CSns,CSs] = obscalnewtest2015( betasigma3, data , basellh, demandcal, p0, rounderr,WFcal)

% [mugamma0 alpha-1 beta gammaishape gammaimean eta-1 r lamda1 lamda2 betacond betapop betalocal olp delta c]


gamma0 = betasigma3(ones(data.N,1))' .* (data.numlist.^betasigma3(8)./mean(data.numlist.^betasigma3(8))); %(data.cdindex,1)
alpha = repmat((betasigma3(2)+1)*betasigma3(14)^betasigma3(15),data.N,1);
beta = betasigma3(3)./betasigma3(14)^betasigma3(15);
m = betasigma3(4);
gammascale =  betasigma3(5)./m.*(data.numlist.^ betasigma3(9)./mean(data.numlist.^ betasigma3(9))).*exp(betasigma3(12).*data.localint);
eta = betasigma3(6)+1; % betasigma3(6*ones(data.N,1))'+1+betasigma3(11).*data.popular;
r = betasigma3(7);
betacond = betasigma3(10);
olp = betasigma3(13);
delta = betasigma3(14);
naturaldisappear = betasigma3(16);

%%%%%%%%%%
%Solve for demand and its first and second order derivatives
%%%%%%%%%%

[D0,dD0,d2D0,sumpind,expU]=demandshopper(alpha,beta,p0- betacond.*data.conditiondif./alpha,data); % 


%% calculate gamma1
p = data.p - rounderr;
Dm=delta.*(p.^(-eta));
dDm=delta.*(-eta).*(p.^(-eta-1));
d2Dm=delta.*(eta).*(eta+1).*(p.^(-eta-2));

%Solve for the gamma that rationalize the price choice
[gamma1,l1]=solvegammaPar(Dm,D0-rounderr.*dD0+0.5.*d2D0.*rounderr^2,dDm,dD0-rounderr.*d2D0,gamma0,p,r);


%% calculate gamma2

p = data.p + rounderr;
Dm=delta.*(p.^(-eta));
dDm=delta.*(-eta).*(p.^(-eta-1));
d2Dm=delta.*(eta).*(eta+1).*(p.^(-eta-2));

%Solve for the gamma that rationalize the price choice
[gamma2,l2]=solvegammaPar(Dm,D0+rounderr .* dD0+0.5.*d2D0 .* rounderr^2,dDm,dD0+rounderr.*d2D0,gamma0,p,r);

SOC = r.*p.*(gamma2.*d2Dm + gamma0.*d2D0) + 2*(r+gamma2.*Dm + gamma0.*(D0+rounderr.*dD0+0.5.*d2D0.*rounderr^2)).*(gamma2.*dDm + gamma0.*(dD0+rounderr.*d2D0));








% pHyp = (data.p - p0).*2;
% DmHyp=delta.*(pHyp.^(-eta));
% dDmHyp=delta.*(-eta).*(pHyp.^(-eta-1));
% 
% gamma4=solvegammaPar(DmHyp,0,dDmHyp,0,0,pHyp,r);




%% check for validity

    
%     palt = [0.7:0.1:1.5 1.7:0.2:5 6:15];
%     lp = length(palt);
%     expUalt = repmat(expU,1,lp) .* min(exp((repmat(data.p,1,lp)+repmat(-palt,data.N,1)).*alpha(1)),1e5);
% %   expUalt = repmat(expU.*exp(data.p.*alpha(1)),1,lp).*repmat(exp(palt.*-alpha(1)),data.N,1);
% %   changepexp = [zeros(data.N,1) repmat(exp(palt.*-alpha(1)),data.N,1) - repmat(exp(data.p.*-alpha),1,lp)];
%     S = repmat(gamma1,1,lp+1).*[Dm delta.*exp(-eta*log(palt))] + repmat(gamma0,1,lp+1).* [expU expUalt] ./ [sumpind repmat(sumpind - expU,1,lp) + expUalt];
%     profittest = [data.p repmat(palt,data.N,1)].*(S./(r + S));
%     check = profittest(:,1) >= (max(profittest,[],2)-1e-12);
%     gamma1(~check,1) = -2;
   

% gammadisplayindex = ((real(gamma1)<0 & real(gamma2)>=0) | (imag(gamma1)~=0) | real(gamma2)<real(gamma1) | SOC>0);
% % disp([data.p(gammadisplayindex,1) gamma1(gammadisplayindex,1) gamma2(gammadisplayindex,1)-gamma1(gammadisplayindex,1) SOC(gammadisplayindex,1)]);
% %     gamma1(dgammaidg<0,1) = -3;    
% %     gamma1(SOC>0,1) = -2; 
% %     gamma1(~l1,1) = -1;
% 
% save diagnewtest.mat
% 
% 
% dsfds

gamma1((imag(gamma1)~=0),1) = real(gamma1((imag(gamma1)~=0),1));
gamma2((imag(gamma2)~=0),1) = 0;
gamma2(SOC>0,1) =0;
gamma2 = max(gamma2,0);

profit2 = p./(r./(gamma2.*Dm+gamma0.*D0)+1);
profitH = ((r*(eta-1)/delta)^(-1/eta)/(1/(eta-1)+1)).*gamma2.^(1/eta);
gamma3=solvegammaPar(Dm,0,dDm,0,0,p,r);
gamma2(profitH>profit2,1) = gamma3(profitH>profit2,1);

gamma1 = min(max(gamma1,0),gamma2);


%% demands functions for true price

% p = data.p;
% Dm=delta.*(p.^(-eta));
% dDm=(-eta).*(p.^(-eta-1));
% d2Dm=(eta).*(eta+1).*(p.^(-eta-2));
% 
% [D0,dD0,d2D0,sumpind,expU]=demandshopper(alpha,beta,p0- betacond.*data.conditiondif./alpha,data); % 
% 

% dDm = delta.*dDm;
% d2Dm = delta.*d2Dm;


 % SD = gamma0.*D0;
% NSD = gamma2.*Dm;
% NSD(gamma2<0,:) = gammascale(gamma2<0,:).*m.*Dm(gamma2<0,:);
% NSDr = gammascale.*m.*Dm;

% %calculate dg/dgammai
% dgdgammai = caldgdgammai(Dm,D0,dDm,dD0,gamma0,gamma1,l1,d2Dm,d2D0,data.p,r);
% dgammaidg= ones(data.N,1)./dgdgammai;






% if length(gamma2)~=length(gammascale) || length(gamma1)~=length(gammascale) || length(m)~=1

% end
Dm=delta.*(data.p.^(-eta));

%%
if demandcal == 1

% D0 =demandshopper(alpha,beta,p0- betacond.*data.conditiondif./alpha,data); % 

    
    demandlh = (data.disappear>0)-(2.*data.disappear -1).*exp(-0.166666667.*(gamma0.*D0 + (gamma2+gamma1)./2.*Dm)).*naturaldisappear;
    demandlhol = (data.disappear>0)-(2.*data.disappear - 1 ) ... Next line starts probability of nondisappear
        .*exp(-0.166666667.*(gamma0.*D0)) ...due to shopper demand. Next line is due to nonshopper demand (taken expectation wrt to gammai)
        .*(1 + gammascale * 0.166666667.* Dm).^-m.*naturaldisappear;
    
%     lip_o = interp1(pdfdata1,pdfdata2,gamma1./gammascale)./gammascale.*dgammaidg.*demandlh.^3; % ;
%     lip_o = exp(-gamma1./gammascale)./gammascale.*dgammaidg.*demandlh.^3; % ;
    lip_o = min((gamcdf(gamma2,m,gammascale) - gamcdf(gamma1,m,gammascale)),1).*demandlh.^3; % ;
    lip_ol = basellh.*  demandlhol.^3; %;  % price likelihood. Next line starts disappear likelihood
else
%    lip_o = interp1(pdfdata1,pdfdata2,gamma1./gammascale)./gammascale.*dgammaidg; % ;
%     lip_o = exp(-gamma1./gammascale)./gammascale.*dgammaidg; % ;
    lip_o = min(gamcdf(gamma2, m,gammascale) - gamcdf(gamma1,m,gammascale),1);
    lip_ol = basellh; %;  % price likelihood. Next line starts disappear likelihood
end

% lip_o(gamma1<0,1) = 0;
% lip_o(lip_o<0,1) = 0;
liptemp = (1-olp).*lip_o + olp.*lip_ol;
olppost = olp.*lip_ol./liptemp;
lip = log(liptemp);

pi = zeros(data.N,1);
CSns = zeros(data.N,1);
CSs = zeros(data.N,1);
if WFcal==1
    [pi,CSns,CSs] = welfaresimple(gamma1,gamma2,gammascale.*m,gamma0,olppost,Dm,D0,p0,data,[alpha(1) beta eta r]);
end

% if isreal(lip)==0
%     save diagwrong.mat
%     dsfds
% end

% if lip(262)>-2
%     save diag.mat
% end

%  save diagnew1.mat

% if length(data.p) == length(data.cdindex)
%     save('July25check')
% end
% if demandcal == 1
%     save('July23check2')
% end

% save('Aug19check')

% 
% if length(data.p)>10000 && betasigma3(1)<4 && betasigma3(14)<3
%     save('Sep8diag')
% end

end

