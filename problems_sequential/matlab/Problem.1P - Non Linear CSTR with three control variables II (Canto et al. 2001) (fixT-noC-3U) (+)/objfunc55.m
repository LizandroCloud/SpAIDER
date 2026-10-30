 function [ J, dJ ] = objfunc( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE )
 if nargout == 1
 f = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE);
 J = -f(end,8);
 else
 [f,df] = fun( x0, ns, ts, ws, ne,n_s,n_c, w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE);
 J = -f(end,8);
 dJ = -df(end,8,:);
 end
 end