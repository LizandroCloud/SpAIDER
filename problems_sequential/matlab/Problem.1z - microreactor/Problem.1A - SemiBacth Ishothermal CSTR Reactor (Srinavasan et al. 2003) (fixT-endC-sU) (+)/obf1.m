function [J] = obf1(f,t_ult,t_out,y,stages)
k1 = 0.053     ;
k2 = 0.128     ;
cb_in = 5      ;
Vo = 1         ;
lamb=1^-8      ;
cao=0.72       ;
cbo=0.05       ;
cd_mais = 0.15 ;
cb_mais = 0.025;

J=-cao*Vo + f(end,1)*f(end,3) ;
end