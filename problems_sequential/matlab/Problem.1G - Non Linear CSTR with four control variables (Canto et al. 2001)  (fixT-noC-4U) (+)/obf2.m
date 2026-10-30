function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)

J=(-f(end,8));
dJ=+(0)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+(0)*df(end,6,:)+(0)*df(end,7,:)+(-1)*df(end,8,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws))+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end