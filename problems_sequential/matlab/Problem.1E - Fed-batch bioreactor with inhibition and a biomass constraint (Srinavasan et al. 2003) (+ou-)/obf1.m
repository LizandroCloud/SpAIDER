function [J] = obf1(f,t_ult,t_out,y,stages)
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
end