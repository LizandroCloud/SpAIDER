function [u] = u_reg(t,ts,w,ks)
lambd=0.01;
un = [];
%for t=1:250
   % w = [ 0.5 -0.4 0.3 0.1 -0.2 0.8 -0.5 0.6];
for i=1:length(w)
    if i==1
     un(1)=0;
     h = 0.5*(1 + tanh( (t+ts(i+1))/lambd ) );  % (função de regularizacao)...
%      g = 0.5*(1 + ( (t+ts(i+1))/sqrt((t+ts(i+1))^2 + lambd^2)) );
    else
     h = 0.5*(1 + tanh( (t-ts(i) )/lambd));  % (função de regularizacao)...
%      g = 0.5*(1 + ( (t-ts(i+1))/sqrt((t-ts(i+1))^2 + lambd^2)) );
    end
    if i==1
        u = w(i)*h;
    else
        u = (w(i)-w(i-1))*h + un(i-1);
    end
    
    un(i) = u;
end
%end

