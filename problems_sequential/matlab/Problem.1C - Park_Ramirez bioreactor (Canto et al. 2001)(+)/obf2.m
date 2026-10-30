function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages,ne)

J=-f(end,1)*f(end,5);
dJ=+(-f(end,5))*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(-f(end,1))*df(end,5,:)+0*ones(1,1,length(ws));
end