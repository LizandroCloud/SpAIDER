function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
k1 = 0.053     ;
k2 = 0.128     ;
cb_in = 5      ;
Vo = 1         ;
lamb=1^-8      ;
cao=0.72       ;
cbo=0.05       ;
cd_mais = 0.15 ;
cb_mais = 0.025;

J=f(end,2) ;
dJ=+(0)*df(end,1,:)+(1)*df(end,2,:)+(0)*df(end,3,:)+0*ones(1,1,length(ws));
end