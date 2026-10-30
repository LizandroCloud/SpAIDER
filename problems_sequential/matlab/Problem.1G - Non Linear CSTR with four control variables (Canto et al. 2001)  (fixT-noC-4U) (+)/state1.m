function [res,ires] = state1(t,x,xp,rpar)
global DEBUG;
global est;
% Parâmetros do modelo... (ver artigo do Schlegel (2004).)
rpar = rpar';
 w(:,1)=rpar(1:((length(rpar)-3)-1)/2,1); 
  ks =  rpar( length( ws(:,1) ) + 1); 
 ts =  rpar(((length(rpar)-2)-2)/2 + 2:length(rpar)-2,1);
 reg = rpar(length(rpar)-1);
 tal = rpar(length(rpar)-0);

% Condição para a regularização...
if reg==1
    
u = u_reg1(t,ts,w,ks);
else
    if est==0
        u = u_reg0(t,ts,w,ks);
    else
        u=w(ks);
    end
end


y(1)=((u(1) +u(2)+u(4)));

dx(1)=(u(4) - y(1)*x(1) - 17.6*x(1)*x(2) - 23*x(1)*x(6)*u(3)                                        )*tal;
dx(2)=(u(1) - y(1)*x(2) - 17.6*x(1)*x(2) - 146*x(2)*x(3)                                            )*tal;
dx(3)=(u(2) - y(1)*x(3)-73*x(2)*x(3)                                                                )*tal;
dx(4)=(-y(1)*x(4)+35.2*x(1)*x(2)-51.3*x(4)*x(5)                                                     )*tal;
dx(5)=(-y(1)*x(5)+219*x(2)*x(3)-51.3*x(4)*x(5)                                                      )*tal;
dx(6)=(-y(1)*x(6)+102.6*x(4)*x(5)-23*x(1)*x(6)*u(3)                                                 )*tal;
dx(7)=(-y(1)*x(7)+46*x(1)*x(6)*u(3)                                                                 )*tal;
dx(8)=(5.8*(y(1)*x(1)-u(4))-3.7*u(1)-4.1*u(2)+y(1)*(23*x(4)+11*x(5)+28*x(6)+35*x(7))-5*u(3)^2 -0.099)*tal;

% res(1)= ( -k1*x(1)*x(2)-( u/x(3) )*x(1)) - xp(1);
% res(2)= ( -k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u/x(3))*(cb_in - x(2)))) - xp(2); 
% res(3)= u - xp(3);
res=dx;

ires = 0;

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
