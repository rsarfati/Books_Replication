function [pi,CSns,CSs] = welfaresimple(gamma1,gamma2,gammascale,gamma0,olppost,Dm,D0,pdif,data,scalarparas)
% scalarpars = [beta,eta,r,delta]

alpha = scalarparas(1);
beta = scalarparas(2);
eta = scalarparas(3);
r = scalarparas(4);


gamma1ave = 0.5*(gamma1+gamma2);

pi_o  = data.p.*(gamma1ave .*Dm + gamma0.*D0)./(r + gamma1ave .*Dm + gamma0.*D0);
pi_ol = data.p.*(gammascale.*Dm + gamma0.*D0)./(r + gammascale.*Dm + gamma0.*D0);
pi = olppost.*pi_ol + (1-olppost).*pi_o;

CSns_o  = (1/(eta-1)).*gamma1ave .*Dm.*data.p./(r + gamma1ave.*Dm + gamma0.*D0);
CSns_ol = (1/(eta-1)).*gammascale.*Dm.*data.p./(r + gammascale.*Dm + gamma0.*D0);
CSns = olppost.*CSns_ol + (1-olppost).*CSns_o;

CSgain = zeros(data.N,1);
N=10000;

for k = 1:data.M

    mktsize = data.cdindex(k)-data.first(k)+1;
    
    randomprice = repmat(-[pdif(data.first(k):data.cdindex(k),1) ; -beta],1,N) ...
        - evrnd(0,1,mktsize+1,N)./alpha(1); 
    
    [best,bestindex] = max(randomprice);
    
    temp = sparse([bestindex mktsize+1], 1:N+1, [best - randomprice(end,:) 1]);
    
    CSgain(data.first(k):data.cdindex(k),1) = sum(temp(1:mktsize,1:N),2) ./ (sum(temp(1:mktsize,1:N)>0,2)+1e-5);

end

CSs_o  = gamma0 .* D0 ./ (r + gamma1ave  .* Dm + gamma0 .* D0) .* CSgain;
CSs_ol = gamma0 .* D0 ./ (r + gammascale .* Dm + gamma0 .* D0) .* CSgain;
CSs = olppost.*CSs_ol + (1-olppost).*CSs_o; 

end

