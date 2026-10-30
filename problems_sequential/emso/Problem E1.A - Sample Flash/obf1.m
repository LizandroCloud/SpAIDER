function [J] =  obf1(f,t_ult,t_out,y)

 J=f(end,1) + 0*t_ult/5000;
%  J=1*abs(f(end,1)-440) + 0*t_ult/5000;
end