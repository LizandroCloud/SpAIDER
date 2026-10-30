function [ dx ] = state( t, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points)

global est;

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

end
