function [J] = obf1(f,t_ult,t_out,y,stages)
k1 = 1;

J=- f(end,1) ;
end