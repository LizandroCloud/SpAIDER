function [J] = obf1(f,t_ult,t_out,y,stages)
Fout = 0.0   ;
Eenz = 1.937 ;
Csol = 80;   ;
Gsol = 0     ;
k = 0.0257458;
Km = 11.6035;;
Kig = 5.38839;

 J=-f(end,2)*f(end,4) + t_ult;
% J=-(100*f(end,2))*f(end,4)/(1.111*(80*0.1 + y(end,5)*y(end,1).*t_ult));

end