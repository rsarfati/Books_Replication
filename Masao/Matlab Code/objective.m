function f = objective(x,x0,distpara0,gamma0vec,deltavec,data12,data09,bp)
xx = [x(1:2) x0(3) x(3:5) x0(7) x(6:(length(x0)-2))];
if xx(4) < 0 || xx(5)<0 || xx(14) >1 || xx(14) < 0
    f = Inf;
else
    f = fullmodelllhWFAug22newtest2015(xx,distpara0,gamma0vec,deltavec,data12,data09,bp);
end
end