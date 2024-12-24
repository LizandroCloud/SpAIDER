function [u] = u_reg0(t,ts,w,ks,ns)
 lambd=0.01;
un = [];
% perfil de u
% if t==125
%    pause
% end
for i=1:length(w)
    if i==1
     un(1)=0;
     h = heaviside(t+ts(i+1));  % (função de regularizacao)...
    %  g = 0.5*(1 + ( (t+ts(i+1))/sqrt((t+ts(i+1))^2 + lambd^2)) );
    else
     h = heaviside( (t-ts(i) ));  % (função de regularizacao)...
      %g = 0.5*(1 + ( (t-ts(i+1))/sqrt((t-ts(i+1))^2 + lambd^2)) );
    end
    if i==1
        %u = w(i)*g;
        u = w(i)*h;
    else
        u = (w(i)-w(i-1))*h + un(i-1);
     %   u = (w(i)-w(i-1))*g + un(i-1);
    end
    
    un(i) = u;
end
%end

