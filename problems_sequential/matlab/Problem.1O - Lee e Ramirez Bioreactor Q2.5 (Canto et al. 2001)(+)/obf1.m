function [J] = obf1(f,t_ult,t_out,y,stages)

J=(-f(end,1)*f(end,4) + 0*f(end,8));
end