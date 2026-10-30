function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
k1 = 0.8   ;
k2 = 0.02  ;
k3 = 0.003 ;
N = 1      ;
R = 8.314  ;
T = 400    ;
V = 1      ;
Pm = 340000;

J=-f(end,3) ;
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(-1)*df(end,3,:)+0*ones(1,1,length(ws));
end