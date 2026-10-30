function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
Feed =  (24.0)/60.0 ;
x_Feed = 0.5        ;
D=u(2)*Feed         ;
vol=1.6             ;
atray=0.25          ;
acond=0.5           ;
areb=1.0            ;

J=( f(end,1) - 0.99)^2 ;
dJ=+(2.0*f(end,1) - 1.98)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+(0)*df(end,6,:)+(0)*df(end,7,:)+(0)*df(end,8,:)+(0)*df(end,9,:)+(0)*df(end,10,:)+(0)*df(end,11,:)+(0)*df(end,12,:)+(0)*df(end,13,:)+(0)*df(end,14,:)+(0)*df(end,15,:)+(0)*df(end,16,:)+(0)*df(end,17,:)+(0)*df(end,18,:)+(0)*df(end,19,:)+(0)*df(end,20,:)+(0)*df(end,21,:)+(0)*df(end,22,:)+(0)*df(end,23,:)+(0)*df(end,24,:)+(0)*df(end,25,:)+(0)*df(end,26,:)+(0)*df(end,27,:)+(0)*df(end,28,:)+(0)*df(end,29,:)+(0)*df(end,30,:)+(0)*df(end,31,:)+(0)*df(end,32,:)+(0)*df(end,33,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end