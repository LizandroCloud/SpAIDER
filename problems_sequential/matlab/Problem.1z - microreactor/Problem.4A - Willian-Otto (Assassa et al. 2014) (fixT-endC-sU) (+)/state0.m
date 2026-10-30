function [res,ires] = state0(t,x,xp,rpar)
global DEBUG;
global est;
% global te=1;
global u;
global ts
global ks
tal = 1;
for i=1:8
if (x(i)<1e-2)
    x(i) = 0;
end
end   

% Parâmetros do modelo... (ver artigo do Schlegel (2004).)
%  rpar = rpar';
%  ws(:,1)=rpar(1:((length(rpar)-3)-2)/2,1); 
%  ks =  rpar( length( ws(:,1) ) + 1); 
%  ts =  rpar(length( ws(:,1) ) + 1 + 2:length(rpar)-2,1);
%  reg = rpar(length(rpar)-1);
%  tal = rpar(length(rpar)-0);


% [u] = ver_reg_est1(ws,1,t);


E1 = -6666.7  ;
E2 = -8333.3  ;
E3 = -11111   ;
k1 = 1.6599e6 ;
k2 = 7.2117e8 ;
k3 = 2.6745e12;
dh1 = -263.8  ;
dh2 = -158.3  ;
dh3 = -226.3  ;
cp = 4.184    ;
cpx = 5554.1  ;
cex = 125.91  ;
Tin = 35      ;
Aref = 9.2903 ;
mref = 2.1052 ;
U = 2.3082e-1 ;
Ma = 100      ;
Mb = 100      ;
Mc = 200      ;
Me = 200      ;
Mg = 300      ;
Mp = 100      ;
roa = 1       ;
rob = 0       ;
roc = 0       ;
roe = 0       ;
rog = 0       ;
rop = 0       ;
tro = 65      ;
Vo = 2        ;
Fscal = 1000  ;
tal = 1000;
y(1)=(exp(-E1/(x(7) + 273.15)));
y(2)=(exp(-E2/(x(7) + 273.15)));
y(3)=(exp(-E3/(x(7) + 273.15)));
% 
% dx(1)=((-x(1)*u(1)/(Fscal*x(8)))-k1*y(1)*x(1)*x(2)                                                                                                                 )*tal;
% dx(2)=(((1-x(2))*u(1)/(Fscal*x(8)))                                                                                           )*tal;
% dx(3)=((-x(3)*u(1)/(Fscal*x(8))) + (Mc*k1*y(1)*x(1)*x(2)/Ma) - (Mc*k2*y(2)*x(2)*x(3)/Mb) - k3*y(3)*x(3)*x(6)                                                       )*tal;
% dx(4)=((-x(4)*u(1)/(Fscal*x(8))) + (Me*k2*y(2)*x(2)*x(3)/Mb)                                                                                                       )*tal;
% dx(5)=((-x(5)*u(1)/(Fscal*x(8))) + (Mg*k3*y(3)*x(3)*x(6)/Mc)                                                                                                       )*tal;
% dx(6)=((-x(6)*u(1)/(Fscal*x(8))) + k2*y(2)*x(2)*x(3) - (Mp*k3*y(3)*x(3)*x(6)/Mc)                                                                                   )*tal;
% dx(7)=(((u(1)*(Tin - x(7)))/(Fscal*x(8))) - (dh1*k1*y(1)*x(1)*x(2)/cp) - (dh2*k2*y(2)*x(2)*x(3)/cp) - (dh3*k3*y(3)*x(3)*x(6)/cp) - (U*Aref/(cp*mref))*(x(7) - 1000*u(2)))*tal;
% % dx(7)= 0;
% dx(8)=(u(1)/Fscal                                                                                                                                                  )*tal;
% % dx(9)=(-cpx*(k2*y(2)*x(2)*x(3) - 0.5*k3*y(3)*x(3)*x(6))*x(8) - cex*(2*k2*y(2)*x(2)*x(3)*x(8))     )*tal;
% dx=transp(dx);


k1 = 1659900.0;
k2 = 721170000.0;
k3 = 1442340000.0;
k4 = 1337250000000.0;
k5 = 4011750000000.0;
k6 = 2674500000000.0;
k7 = 3319800.0;
k8 = 104656218.9;
k9 = 27285184270.0;
k10 = 144655676400000.0;
k11 = 181605029400.0;
l1 = 0.0002434546857;
l2 = 1000;
Tin = 35;
b2 = 6.6667;

dx(1)=((-x(1)*u(1)/(Fscal*x(8)))-k1*y(1)*x(1)*x(2)                                                                                                                 )*tal;
dx(2)=(((1-x(2))*u(1)/(Fscal*x(8)))-k1*y(1)*x(1)*x(2)- k2*y(2)*x(2)*x(3)                                                                                           )*tal;
dx(3)=((-x(3)*u(1)/(Fscal*x(8))) + (k7*y(1)*x(1)*x(2)) - (k3*y(2)*x(2)*x(3)) - k6*y(3)*x(3)*x(4))*tal; 
dx(4)=((-x(4)*u(1)/(Fscal*x(8))) + k2*y(2)*x(2)*x(3) - (k4*y(3)*x(3)*x(4)) )*tal;
% dx(4)=((-x(4)*u(1)/(Fscal*x(8))) + (Me*k2*y(2)*x(2)*x(3)/Mb)                                                                                                       )*tal;
dx(5)=((-x(5)*u(1)/(Fscal*x(8))) + (k3*y(2)*x(2)*x(3))                                                                                                       )*tal;
dx(6)=((-x(6)*u(1)/(Fscal*x(8))) + k5*y(3)*x(3)*x(4)                                                                                    )*tal;
dx(7)=(((u(1)*(Tin - x(7)))/(Fscal*x(8))) - (k8*y(1)*x(1)*x(2)) - (k9*y(2)*x(2)*x(3)) - (k10*y(3)*x(3)*x(4)) - l1*x(7) + l2*u(2))*tal;
dx(8)=(u(1)/Fscal     )*tal   ;                


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
