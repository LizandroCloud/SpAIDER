function [res,ires] = state1(t,x,xp,rpar)
global DEBUG;
global est;
% Parâmetros do modelo... (ver artigo do Schlegel (2004).)
rpar = rpar';
 w(:,1)=rpar(1:((length(rpar)-3)-1)/2,1); 
 ks =  rpar(((length(rpar)-2)-2)/2 + 1); 
 ts =  rpar(((length(rpar)-2)-2)/2 + 2:length(rpar)-2,1);
 reg = rpar(length(rpar)-1);
 tal = rpar(length(rpar)-0);

% Condição para a regularização...
if reg==1
    
u = u_reg1(t,ts,w,ks);
else
    if est==0
        u = u_reg0(t,ts,ws,ks);
    else
        u=ws(ks);
    end
end


% [u] = ver_reg(ws,t_out,te,ks,ns,reg,est);
% [uref] = ver_reg(utarg,t_out,tstarg,ks,nstarg,reg,est);
mim = 0.53;
km = 1.2  ;
ki = 22   ;
Yx = 0.4  ;
Yp = 1    ;
v = 0.5   ;
Sin = 20  ;
xo = 1    ;
So = 0    ;
Po = 0    ;
Vo = 2;   ;
xmax= 3   ;
smax= 5   ;

y(1)=(mim*x(2) / (km + x(2)+ ( (x(2)^2)/ki )));

res(1)=(y(1)*x(1) - u(1)*x(1)/x(4)                       )*tal - xp(1);
res(2)=(- y(1)*x(1)/Yx - v*x(1)/Yp + u(1)*(Sin-x(2))/x(4) )*tal - xp(2);
res(3)=( v*x(1) - u(1)*x(3)/x(4)                         )*tal - xp(3);
res(4)=(u(1)                                             )*tal - xp(4);

% res(1)= ( -k1*x(1)*x(2)-( u/x(3) )*x(1)) - xp(1);
% res(2)= ( -k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u/x(3))*(cb_in - x(2)))) - xp(2); 
% res(3)= u - xp(3);


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
