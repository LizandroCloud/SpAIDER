function [ dx ] = state( t, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,points,ns_a,frozen,ts_frozen)

global est;


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

end
