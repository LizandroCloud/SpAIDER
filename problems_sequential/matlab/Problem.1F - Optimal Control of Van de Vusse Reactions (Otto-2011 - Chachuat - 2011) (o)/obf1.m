function [J] = obf1(f,t_ult,t_out,y,stages)
V=3.14     ;
U=642.18   ;
A=25.1     ;
k1=1.287e12;
k2=1.287e12;
k3=9.043e9 ;
E1=9758.3  ;
E2=9758.3  ;
E3=8560    ;
dH1=4200   ;
dH2=-11000 ;
dH3=-41850 ;
ro=1000    ;
cp=2.8121  ;
caf=5.5    ;
cbf=0      ;
ccf=0      ;
cdf=0      ;
tf=378     ;

J=-f(end,2);
end