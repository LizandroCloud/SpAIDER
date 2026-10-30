function [ dx ] = sens( t, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,check,uws)

global est;

u(1:length(ws(:,1)) )=0;

mim = 0.53; % l/(h)
km = 1.2;   % (g/L)
ki = 22; % (g/L)
Yx = 0.4;
Yp = 1;
v = 0.5; %L/h;
Sin = 20; % g/l
xo = 1; %g/L;
So = 0; %g/L;
Po = 0; %g/L;
Vo = 2; %L

% Regularization condition...
for i = 1:length(ws(:,1)) 
if reg==1
     u(i) = u_reg1(t,te(i,:),ws(i,:),ks,ne(i));
else
    if est==0
        u(i) = u_reg0(t,te(i,:),ws(i,:),ks,ne(i));
       else
        u(i) = ws(i);
    end
 end
end

mi=mim*x(2) / (km + x(2)+ ( (x(2)^2)/ki ));


dx(1)= mi*x(1) - u*x(1)/x(4); % X
dx(2)= - mi*x(1)/Yx - v*x(1)/Yp + u*(Sin-x(2))/x(4); % S
dx(3)=  v*x(1) - u*x(3)/x(4);  % P
dx(4) = u; %V
dx=dx';

% 
%  syms x2 x1 x3 x4 u mi Yx v Sin Yp
%  eq = [(mim*x2 / (km + x2+ ( (x2^2)/ki )))*x1 - u*x1/x4, - (mim*x2 / (km + x2+ ( (x2^2)/ki )))*x1/Yx - v*x1/Yp + u*(Sin-x2)/x4,  v*x1 - u*x3/x4, u];
%  [jac,j1] = jacmat( eq , [x1,x2,x3,x4,u], x );


 
 
jc1 = [   53/100*x(2)/(6/5+x(2)+1/22*x(2)^2)-u/x(4),             53/100/(6/5+x(2)+1/22*x(2)^2)*x(1)-53/100*x(2)/(6/5+x(2)+1/22*x(2)^2)^2*x(1)*(1+1/11*x(2)),                                                                                      0,                                                                              u*x(1)/x(4)^2,                                                                                 -x(1)/x(4)
          -53/100*x(2)/(6/5+x(2)+1/22*x(2)^2)/Yx-v/Yp, -53/100/(6/5+x(2)+1/22*x(2)^2)*x(1)/Yx+53/100*x(2)/(6/5+x(2)+1/22*x(2)^2)^2*x(1)/Yx*(1+1/11*x(2))-u/x(4),                                                                                      0,                                                                       -u*(Sin-x(2))/x(4)^2,                                                                            (Sin-x(2))/x(4)
          v,                                                                                      0,                                                                                  -u/x(4),                                                                              u*x(3)/x(4)^2,                                                                                 -x(3)/x(4)
          0,                                                                                      0,                                                                                      0,                                                                                      0,                                                                                      1];
 



jk=double(jc1);
j0=0;
jk(:,1:n_s)=jc1(:,1:n_s);
jk(:,n_s+1:n_s+n_c)=0;



for j=1:n_c
    
%%%%%      tx = [t ts_frozen(1,2:ns_a(j)+1)];
%     tx = [t ts_frozen(check(j,1:length(check)))];
%     tx = sort(tx);
%     posi = find(tx == t,1,'last');
    
    for is = 1:ns_a(j) %length(ws)
%         if ((posi)== is)
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