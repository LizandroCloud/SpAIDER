function [J dJ] = obf2(f,t_ult,t_out,df,ws,y)

J=f(end,3);
dJ=(0)*df(end,1,:)+(0)*df(end,2,:)+(1)*df(end,3,:)+0*ones(1,1,length(ws));
end