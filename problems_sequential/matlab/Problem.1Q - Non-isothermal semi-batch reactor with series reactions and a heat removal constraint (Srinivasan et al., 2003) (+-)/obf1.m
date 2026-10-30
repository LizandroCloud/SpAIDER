function [J] = obf1(f,t_ult,t_out,y,stages)
k10  = 4    ;
k20 = 800   ;
E1 = 6e3    ;
R = 8.31    ;
dH1 = -3e4  ;
dH2 = -1e4  ;
qr_maf=1.5e5;
V_maf=1.1   ;
E2 = 2e4    ;
cb_in = 20  ;

J=f(end,3)*f(end,4);
end