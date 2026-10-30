 function [ J, dJ ] = obj( x0, ns, ts, ws, optODE )
 if nargout == 1
 f = fun( x0, ns, ts, ws, optODE );
 J = f(3);
 else
 [f,df] = fun( x0, ns, ts, ws, optODE );
 J = f(3);
 dJ = df(3,:)';
 end
 end