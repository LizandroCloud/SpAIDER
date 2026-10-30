function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)

J=( -f(end,2) );
dJ=+(0)*df(end,1,:)+(-1)*df(end,2,:)+0*ones(1,1,length(ws));
end