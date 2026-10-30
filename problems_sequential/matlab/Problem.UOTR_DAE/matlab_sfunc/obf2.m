function [J,gradJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
gradJ=[];

%Objective
J = f(end,end);

%Gradient
gradJ = df(end,end,:); gradJ=gradJ(:)';
end
