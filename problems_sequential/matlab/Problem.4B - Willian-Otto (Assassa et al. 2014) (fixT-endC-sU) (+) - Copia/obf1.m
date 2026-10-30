function [J] = obf1(f,t_ult,t_out,y,stages)
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
end