function [J] = obf1(f,t_ult,t_out,y,stages)
R = 8.314      ;
n = 2          ;
Ea = 63.717*1e3;
A = 7.5581     ;

J=-f(end,1) + 1e-4*f(end,2);
end