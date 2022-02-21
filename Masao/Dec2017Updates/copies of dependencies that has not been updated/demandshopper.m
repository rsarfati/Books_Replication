function [f1,f2,f3,sum2,expU] = demandshopper( alpha, beta ,p,data)


expU=exp(p.*(-alpha));

%sum utility for each book
% sum1=zeros(M);
% sum1(1)=sum(expU(1:cdindex(1)))+exp(beta);
% for j=2:M
%     sum1(j)=sum(expU(cdindex(j-1)+1:cdindex(j)))+exp(beta);
% end

sum1 = sum(sparse(1:length(data.cdid),data.cdid,data.obsweight.*expU))' + exp(beta.*alpha(1));
sum2 = sum1(data.cdid);

%denominator for searchers
% denom = 1./(sum1(cdid));
% f=expU.*denom;

f1 = expU ./ sum2;
f2 = -alpha.*expU.*(sum2-expU)./(sum2.^2);
f3 = alpha.^2.*expU.*(sum2-expU).*(sum2-2*expU)./(sum2.^3);
%hello= isnan(f);
f1(isnan(f1))=0;
f2(isnan(f2))=0;
f3(isnan(f3))=0;
end