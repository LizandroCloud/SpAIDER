function [J] = obf1(f,t_ult,t_out,y,stages)

J=trapz(0.278 - f(:,1));
end