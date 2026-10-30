function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
M 	= 0.028  ;
Ta = 293    ;
La	= 230.87 ;
Va	= 22.3e-3;
ro	= 923.9  ;
Tw	= 293    ;
Lw	= 2400   ;
Lr	= 40     ;
Ar	= 0.0182 ;
Aw  = 0.0182;
Civ	= 1.6e-3;
Cpc	= 0.5e-2;
Cr 	= 0.8e-4;
rgo = 0     ;
wgc = 1.36e4;
Ps	= 6276e3 ;
Pr	= 18632e3;
vo = 1/ro   ;
R = 8.314   ;
g = 10      ;

J=-f(end,3);
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(-1)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+0*ones(1,1,length(ws));
end