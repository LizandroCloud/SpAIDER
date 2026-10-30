function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)

%Objective
J = f(end,38);

%Gradient
dJ = df(end,38,:); 
dJ=dJ(:)';
end