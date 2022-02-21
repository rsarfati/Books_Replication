function [gamma1,l1] = solvegammaPar(Dm,D0,dDm,dD0,gamma0,p,r)

A=Dm.*(Dm);
B=(dDm).*(r).*(p) +r.*Dm+2*Dm.*gamma0.*D0;
C=r.*p.*gamma0.*dD0+r.*gamma0.*D0+gamma0.*gamma0.*D0.*D0;

gamma1=((-B+(B.^2-4*A.*C).^(0.5))./(2*A));
l1=(real(gamma1)>0)&(imag(gamma1)==0);

% if mean(l1)~=1
%     index = (l1~=1);
%     [A(index,1) B(index,1) C(index,1) gamma1(index,1)]
% end
% gamma1=l1.*gamma1+(1-l1).* -1;
%(-10-10*abs(B.^2-4*A.*C)-10*abs(real(gamma1)));

end




