 function [ J, dJ ] = objfunc( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE )

 
 if nargout == 1
 f = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE);
 J = f(end,5);
%  f(end,4)
 else
 [f,df] = fun(  x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE );

 J = f(end,5);
 dJ = df(end,5,:);
%  f(end,4)
 end
 end