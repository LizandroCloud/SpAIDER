function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
UA = 1560     ;
V = 9450      ;
Vj = 2000     ;
dh = -69919.56;
Rg = 1.987    ;
roj = 1       ;
ros = 0.915   ;
cpj = 4045.7  ;
cps = 1647.27 ;
T0 = 333      ;
Tj0 = 294     ;
ci0 = 0.981   ;
cm0 = 8.63    ;
cb0 = 1.05    ;

J=-f(end,10) ;
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+(0)*df(end,6,:)+(0)*df(end,7,:)+(0)*df(end,8,:)+(0)*df(end,9,:)+(-1)*df(end,10,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end