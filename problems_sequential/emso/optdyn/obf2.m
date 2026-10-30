function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)

J=f(end,1);
dJ=+(1)*df(end,1,:)+0*ones(1,1,length(ws));
end