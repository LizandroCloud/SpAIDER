function [u] = u_regS(t,ts,w,ks,ns)
lambd=0.01;
un = [];
%for t=1:250
   % w = [ 0.5 -0.4 0.3 0.1 -0.2 0.8 -0.5 0.6];
for i=1:ns
    if i==1
     un(1)=0;
     h = 0.5*(1 + tanh( (t+ts(i+1))/lambd ) );  % (regularization)...
     g = 0.5*(1 + ( (t+ts(i+1))/sqrt((t+ts(i+1))^2 + lambd^2)) );
      
    else
     h = 0.5*(1 + tanh( (t-ts(i) )/lambd));  % (regularization)...
     g = 0.5*(1 + ( (t-ts(i))/sqrt((t-ts(i))^2 + lambd^2)) );
    end
    if i==1
        u = w(i)*g;
    else
        u = (w(i)-w(i-1))*g + un(i-1);
    end
    
    un(i) = u;
%     h
end
%end

