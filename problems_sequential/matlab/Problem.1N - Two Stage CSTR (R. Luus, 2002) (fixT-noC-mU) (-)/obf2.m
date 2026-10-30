function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
tau=0.1;

J=f(end,5);
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(1)*df(end,5,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end