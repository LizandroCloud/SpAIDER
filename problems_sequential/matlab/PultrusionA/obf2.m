function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
R = 8.314      ;
n = 2          ;
Ea = 63.717*1e3;
A = 7.5581     ;

J=-f(end,1) + 1e-4*f(end,2);
dJ=+(-1.0)*df(end,1,:)+(0.0001)*df(end,2,:)+0*ones(1,1,length(ws));
end