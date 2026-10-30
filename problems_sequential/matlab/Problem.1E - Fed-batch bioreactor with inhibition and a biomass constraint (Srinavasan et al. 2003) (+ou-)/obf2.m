function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
mim = 0.53;
km = 1.2  ;
ki = 22   ;
Yx = 0.4  ;
Yp = 1    ;
v = 0.5   ;
Sin = 20  ;
xo = 1    ;
So = 0    ;
Po = 0    ;
Vo = 2;   ;
xmax= 3   ;
smax= 5   ;

J=-f(end,3)*f(end,4);
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(-f(end,4))*df(end,3,:)+(-f(end,3))*df(end,4,:)+0*ones(1,1,length(ws));
end