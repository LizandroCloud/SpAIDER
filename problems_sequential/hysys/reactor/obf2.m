function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
k1 = 1;

J=- f(end,1) ;
dJ=+(-1)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+0*ones(1,1,length(ws));
end