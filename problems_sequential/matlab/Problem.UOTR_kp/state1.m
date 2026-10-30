function [res,ires] = state1(t,x,xp,rpar)
global DEBUG;
global est;
% Parâmetros do modelo... (ver artigo do Schlegel (2004).)
 rpar = rpar';
 ws=rpar.ws; 
 ks =  rpar.ks; 
 ts =  rpar.ts;
 reg = rpar.reg;
 tal = rpar.tal;
 n_c =  rpar.nc;

% % Condição para a regularização...
% if reg==1
% 
% u = u_reg1(t,ts,w,ks);
% else
%     if est==0
%         u = u_reg0(t,ts,ws,ks);
%     else
%         u=ws(ks);
%     end
% end


global modelo_u0 modelo_par
global idu

idu=[1;16;17;18];
idu=[1;16; 17];
for j=1:n_c
    ws_reg(j)=ver_reg_est1(ws(j), [],t);
    % ws_reg(j)=u_reg1(t,ts,ws(1:j),ks);
    % ws_reg(j) = u_reg1(t_out,te(j,:),ws(j),ks,ns_a(j));
    % ws_reg=ver_reg_est1(ws(j),yi,t_out);
    % ws_reg(j) = interp1(te(j,:), [ws(j,:),ws(j,end)], t_out, 'cubic');
    % wsb_reg(j) = u_reg1(t_out,te(j,:),ws(j,:),ks,ne(j));
    % det=ws_reg(j) - wsb_reg(j);
    % if det>0.01
    %     det
    % end

end
%Input interpolation
uf=modelo_u0';
for j=1:n_c
    uf(idu(j))=ws_reg(j);
end
uf(idu(1))=uf(idu(1))*1e3;

%States' derivative
dx=uotr_ctrl_Kp(t,x(1:38),uf,1,[],modelo_par);
res = real(dx)*tal - xp;

% res(1)=(y(1)*x(1) - u(1)*x(1)/x(4)                       )*tal - xp(1);
% res(2)=(- y(1)*x(1)/Yx - v*x(1)/Yp + u(1)*(Sin-x(2))/x(4) )*tal - xp(2);
% res(3)=( v*x(1) - u(1)*x(3)/x(4)                         )*tal - xp(3);
% res(4)=(u(1)                                             )*tal - xp(4);

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
