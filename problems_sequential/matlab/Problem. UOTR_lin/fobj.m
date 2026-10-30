function [J,gradJ] = fobj(f,t_ult,t_out,df,ws,y,stages)
%Objective
J = f(end,end);

%Gradient
gradJ = df(end,end,:); gradJ=gradJ(:)';
end