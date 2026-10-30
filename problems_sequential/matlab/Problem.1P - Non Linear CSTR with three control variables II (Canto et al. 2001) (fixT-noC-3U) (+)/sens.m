function [ dx ] = sens( t, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws)

global est;

% dx=zeros(1,length(x));
% if est==0
% koz=1;
% for i=1:length(ne)
%     for j=1:ne(i)+1
%     te(i,j)=ts(koz);
%     koz=koz+1;
%     end
% end
% end
% Regularization condition...
i = 1:length(ws(:,1)) ;
if reg==1
     u(i) = u_reg1(t,te(i,:),ws(i,:),ks,ne(i));
else
    if est==0
        u(i) = u_reg0(t,te(i,:),ws(i,:),ks,ne(i));
       else
        u(i) = ws(i);
    end
 end



q = u(1) +u(2)+u(4);
dx(1)=u(4) - q*x(1) - 17.6*x(1)*x(2) - 23*x(1)*x(6)*u(3);
dx(2)=u(1) - q*x(2) - 17.6*x(1)*x(2) - 146*x(2)*x(3);
dx(3)=u(2) - q*x(3)-73*x(2)*x(3);
dx(4)=-q*x(4)+35.2*x(1)*x(2)-51.3*x(4)*x(5);
dx(5)=-q*x(5)+219*x(2)*x(3)-51.3*x(4)*x(5);
dx(6)=-q*x(6)+102.6*x(4)*x(5)-23*x(1)*x(6)*u(3);
dx(7)=-q*x(7)+46*x(1)*x(6)*u(3);
dx(8)=5.8*(q*x(1)-u(4))-3.7*u(1)-4.1*u(2)+q*(23*x(4)+11*x(5)+28*x(6)+35*x(7))-5*u(3)^2 -0.099;

dx=dx';

% syms x2 x1 x3 x4 x5 x6 x7 x8 u1 u2 u3 u4 q
% eq = [u4 - (u1+u2+u4)*x1 - 17.6*x1*x2 - 23*x1*x6*u3
%     u1 - (u1+u2+u4)*x2 - 17.6*x1*x2 - 146*x2*x3
%       u2 - (u1+u2+u4)*x3-73*x2*x3 
%       -(u1+u2+u4)*x4+35.2*x1*x2-51.3*x4*x5
%     -(u1+u2+u4)*x5+219*x2*x3-51.3*x4*x5
%     -(u1+u2+u4)*x6+102.6*x4*x5-23*x1*x6*u3
%     -(u1+u2+u4)*x7+46*x1*x6*u3
%     5.8*((u1+u2+u4)*x1-u4)-3.7*u1-4.1*u2+(u1+u2+u4)*(23*x4+11*x5+28*x6+35*x7)-5*u3^2 - 0.099];
%    
% [jac,j1] = jacmat( eq , [x1,x2,x3,x4, x5,x6,x7,x8,u1,u2,u3,u4], x );

jc1 = [           -u(1)-u(2)-u(4)-88/5*x(2)-23*x(6)*u(3),                              -88/5*x(1),                                     0,                                     0,                                     0,                             -23*x(1)*u(3),                                     0,                                     0,                                   -x(1),                                   -x(1),                             -23*x(1)*x(6),                                  1-x(1)
                              -88/5*x(2),              -u(1)-u(2)-u(4)-88/5*x(1)-146*x(3),                               -146*x(2),                                     0,                                     0,                                     0,                                     0,                                     0,                                  1-x(2),                                   -x(2),                                     0,                                   -x(2)
                                     0,                                -73*x(3),                       -u(1)-u(2)-u(4)-73*x(2),                                     0,                                     0,                                     0,                                     0,                                     0,                                   -x(3),                                  1-x(3),                                     0,                                   -x(3)
                              176/5*x(2),                              176/5*x(1),                                     0,                   -u(1)-u(2)-u(4)-513/10*x(5),                            -513/10*x(4),                                     0,                                     0,                                     0,                                   -x(4),                                   -x(4),                                     0,                                   -x(4)
                                     0,                                219*x(3),                                219*x(2),                            -513/10*x(5),                   -u(1)-u(2)-u(4)-513/10*x(4),                                     0,                                     0,                                     0,                                   -x(5),                                   -x(5),                                     0,                                   -x(5)
                             -23*x(6)*u(3),                                     0,                                     0,                              513/5*x(5),                              513/5*x(4),                    -u(1)-u(2)-u(4)-23*x(1)*u(3),                                     0,                                     0,                                   -x(6),                                   -x(6),                             -23*x(1)*x(6),                                   -x(6)
                              46*x(6)*u(3),                                     0,                                     0,                                     0,                                     0,                              46*x(1)*u(3),                             -u(1)-u(2)-u(4),                                     0,                                   -x(7),                                   -x(7),                              46*x(1)*x(6),                                   -x(7)
               29/5*u(1)+29/5*u(2)+29/5*u(4),                                     0,                                     0,                     23*u(1)+23*u(2)+23*u(4),                     11*u(1)+11*u(2)+11*u(4),                     28*u(1)+28*u(2)+28*u(4),                     35*u(1)+35*u(2)+35*u(4),                                     0, 29/5*x(1)-37/10+23*x(4)+11*x(5)+28*x(6)+35*x(7), 29/5*x(1)-41/10+23*x(4)+11*x(5)+28*x(6)+35*x(7),                                -10*u(3),  29/5*x(1)-29/5+23*x(4)+11*x(5)+28*x(6)+35*x(7)];
jk=double(jc1);
j0=0;
jk(:,1:n_s)=jc1(:,1:n_s);
jk(:,n_s+1:n_s+n_c)=0;
for j=1:n_c
    
%  %%%    tx = [t ts_frozen(1,2:ns_a(j)+1)];
%     tx = [t ts_frozen(j,(points(j,1:ne(j)-1)))];
%     tx = sort(tx);
%     check = find(tx == t,1,'last');
    
    for is = 1:ns_a(j) %length(ws)
%         if ((check)== is)
%              uws = 1;
%         else   % do contrario...
%              uws = 0;
%         end
                jk(:,n_s+j)=uws(j,is)*jc1(:,n_s+j);
                k=1:(n_s);
                dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0)+ jk(k,n_s+j);
    end 
j0=length(dx)-n_s;
end






% for is = 1:ne
%    
% if is == ks  % ks e o numero de estagios... se for estagio atual tem sensibilidade.
% uws = 1;
% else   % do contrario...
% uws = 0;
% end
% 
% 
% dx(3*is+1)=  (-k1*x(2) - ( u/x(3) ))*x(3*is+1) - k1*x(1)*x(3*is+2) +u*x(1)*x(3*is+3)/(x(3)^2) - uws*x(1)/x(3);
% dx(3*is+2)=  (-k1*x(2) )*x(3*is+1) + ((-u/x(3) )-4*k2*x(2) - k1*x(1))*x(3*is+2) -  u/(x(3)^2)*((cb_in - x(2)))*x(3*is+3) + ((cb_in - x(2))/x(3))*uws;
% dx(3*is+3)=   uws;
% end     
% 
% 
% 
% toc


end
