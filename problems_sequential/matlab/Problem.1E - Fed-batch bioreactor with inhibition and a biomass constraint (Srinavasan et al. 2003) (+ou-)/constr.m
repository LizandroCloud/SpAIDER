function [ c, ceq, dc, dceq ] = constr( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE  )
 
 xmax= 3; %g/L;
 smax=5;
 if nargout == 2   % caso nao tenha analise de sensibilidade...
 [f] = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE );
 ceq = []; 
 for i=1:length(f(:,1));
   c(i) = f(i,1) - xmax;
 end
  k=i;
 for j=1:length(f(:,1));
   c(k) = f(j,2) - smax;
   k=k+1;
 end
 
 
 c=c';
 else              % case of analytic jacobian...
 
 [f,df] = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE );

%   syms f3 f2 f1 cb_in cao cb_in cbo_in Vo
%   eq = [(0.5/f3)*((f1+cb_in-f2)*f3 - (cao+cb_in-cbo)*Vo)];
%   [jac,j1] = jacmat( eq , [f1,f2,f3], f );
 ceq = [];  
 dceq = [];  
 
 for i=1:length(f(:,1));
   c(i) = f(i,1) - xmax;
   dc(:,i) = df(i,1,:);
 end
 k=i;
 for j=k+1:k+length(f(:,1));
   c(j) = f(j,2) - xmax;
   dc(:,j) = df(j,2,:);
 end
 
 ceq = [];  
 dceq = [];  
   
end
end