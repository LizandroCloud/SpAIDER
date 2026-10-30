function [res,ires] = state1(t,x,xp,rpar)
global DEBUG;
global est;
% global te=1;
global u;
global ts
global ks
tal = 1;
% Parâmetros do modelo... (ver artigo do Schlegel (2004).)
%  rpar = rpar';
%  ws(:,1)=rpar(1:((length(rpar)-3)-2)/2,1); 
%  ks =  rpar( length( ws(:,1) ) + 1); 
%  ts =  rpar(length( ws(:,1) ) + 1 + 2:length(rpar)-2,1);
%  reg = rpar(length(rpar)-1);
%  tal = rpar(length(rpar)-0);


% [u] = ver_reg_est1(ws,1,t);


UA = 1560     ;
V = 9450      ;
Vj = 2000     ;
dh = -69919.56;
Rg = 1.987    ;
roj = 1       ;
ros = 0.915   ;
cpj = 4045.7  ;
cps = 1647.27 ;
T0 = 333      ;
Tj0 = 294     ;
ci0 = 0.981   ;
cm0 = 8.63    ;
cb0 = 1.05    ;

y(1)=(9.11e13*exp(-29508/(Rg*x(8)))                                                 );
y(2)=(1.11e5*exp(-27340/(Rg*x(8)))                                                  );
y(3)=(1.01e7*exp(-7067/(Rg*x(8)))                                                   );
y(4)=(2.01e6*exp(-7067/(Rg*x(8)))                                                   );
y(5)=(1.01e7*exp(-7067/(Rg*x(8)))                                                   );
y(6)=(1.01e7*exp(-7067/(Rg*x(8)))                                                   );
y(7)=(6.61e7*exp(-14400/(Rg*x(8)))                                                  );
y(8)=(2.31e9*exp(-18000/(Rg*x(8)))                                                  );
y(9)=(( cm0 - x(2) )/ cm0                                                           );
y(10)=(2.57 - 0.00505*x(8)                                                           );
y(11)=(9.56 - 0.0176*x(8)                                                            );
y(12)=(-3.03 + 0.00785*x(8)                                                          );
y(13)=(1.71e9*exp(-843/x(8) - 2*( y(10)*(y(9)) +  y(11)*(y(9))^2 +  y(12)*(y(9))^3 )));

dx(1)=(( ( u(1)*ci0 - u(2)*x(1) ) / V )  - y(1)*x(1)                                                                                                    )*tal;
dx(2)=(( (u(2) * (cm0 - x(2)) / V )  - y(6)*x(2)*(x(6)+x(7)) )                                                                                          )*tal;
dx(3)=(( (u(2) * (cb0 - x(3)) / V )  - x(3)*( y(4)*x(4) + y(7)*x(6) + y(8)*x(7) ) )                                                                     )*tal;
dx(4)=(( 2*0.57*y(1)*x(1) - x(4)*(y(3)*x(2) + y(4)*x(3)) )                                                                                              )*tal;
dx(5)=(x(3)*(y(4)*x(4) + y(8)*(x(6)+x(7))) - x(5)*(y(5)*x(2) + y(13)*(x(5) + x(6) + x(7) ) )                                                            )*tal;
dx(6)=(2*y(2)*x(2)^3 + y(3)*x(4)*x(2) + x(2)*y(7)*(x(6)+x(7)) + y(6)*x(2)*x(6) - (y(6)*x(2) + y(13)*(x(5) + x(6) + x(7) ) + y(7)*x(2) + y(8)*x(3) )*x(6))*tal;
dx(7)=(y(6)*x(5)*x(2) +  y(6)*x(2)*x(7) - (y(6)*x(2) + y(13)*( x(5) + x(6) + x(7) ) + y(7)*x(2) + y(8)*x(3) )*x(7)                                      )*tal;
dx(8)=(u(2)*(T0 - x(8))/V - (dh*y(6)*x(2)*(x(6)+x(7))/(ros*cps))  - UA*((x(8) - x(9))/(ros*cps*V))                                                      )*tal;
dx(9)=((u(3)*(Tj0 - x(9))/Vj) - UA*(x(8) - x(9))/(roj*cpj*Vj)                                                                                           )*tal;
dx(10)=(y(9)                                                                                                                                             )*tal;
% dx=transp(dx);
res = xp - dx';
ires=0;

if ( DEBUG > 9 )
  buf1 = sprintf('%.4e ',y);
  buf2 = sprintf('%.4e ',yp);
  buf3 = sprintf('%.4e ',res);
  disp(sprintf('(Debug) dydt.m: t=%.2e y=%s yp=%s res=%s', ...
               t,buf1,buf2,buf3));
  if ( DEBUG > 99 )
    cbuf = input(sprintf('\tRETURN to continue : '));
  end
end


end
