function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
Fout = 0.0   ;
Eenz = 1.937 ;
Csol = 80;   ;
Gsol = 0     ;
k = 0.0257458;
Km = 11.6035;;
Kig = 5.38839;

J=-f(end,2)*f(end,4);
dJ=+(0)*df(end,1,:)+(-f(end,4))*df(end,2,:)+(0)*df(end,3,:)+(-f(end,2))*df(end,4,:)+0*ones(1,1,length(ws))+0*ones(1,1,length(ws));
end