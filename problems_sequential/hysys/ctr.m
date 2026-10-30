function [ c, ceq, dc, dceq ] = ctr( x0, ns, ts, ws, optODE )
 if nargout == 2
 f = fun( x0, ns, ts, ws, optODE );
 c = [];
 ceq = [];
 else
 [f,df] = fun( x0, ns, ts, ws, optODE );
 c = f(1:2)';
 dc = df(1:2,:)';
 ceq = [];
dceq = [];
end
end