function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
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
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(f(end,4))*df(end,3,:)+(f(end,3))*df(end,4,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end