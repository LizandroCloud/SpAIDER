function [J, dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)

J=-1+f(end,1)+f(end,2);
dJ=+(1)*df(end,1,:)+(1)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+(0)*df(end,6,:)+(0)*df(end,7,:)+(0)*df(end,8,:)+(0)*df(end,9,:)+(0)*df(end,10,:)+(0)*df(end,11,:)+(0)*df(end,12,:)+(0)*df(end,13,:)+(0)*df(end,14,:)+(0)*df(end,15,:)+(0)*df(end,16,:)+(0)*df(end,17,:)+(0)*df(end,18,:)+(0)*df(end,19,:)+(0)*df(end,20,:)+(0)*df(end,21,:)+(0)*df(end,22,:)+(0)*df(end,23,:)+(0)*df(end,24,:)+(0)*df(end,25,:)+(0)*df(end,26,:)+(0)*df(end,27,:)+(0)*df(end,28,:)+0*ones(1,1,length(ws));
end