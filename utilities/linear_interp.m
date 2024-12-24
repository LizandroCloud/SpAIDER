function [u] = linear_interp(t,ts,w,ks,ns)
% lambd=0.1;
% un = [];
%for t=1:250
   % w = [ 0.5 -0.4 0.3 0.1 -0.2 0.8 -0.5 0.6];
times=find(ts>t);
j=times(1)-1;
sum= w(j); 
% for i=2:2
    coef=(w(j+1)-w(j))/(ts(j+1)-ts(j))*(t-ts(j));
    sum=coef+sum;
% end
   
    u = sum;
end
%end

