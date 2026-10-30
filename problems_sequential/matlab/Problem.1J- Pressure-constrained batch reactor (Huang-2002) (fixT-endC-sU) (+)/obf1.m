function [J] = obf1(f,t_ult,t_out,y,stages)
k1 = 0.8   ;
k2 = 0.02  ;
k3 = 0.003 ;
N = 1      ;
R = 8.314  ;
T = 400    ;
V = 1      ;
Pm = 340000;

J=-f(end,3) ;
end