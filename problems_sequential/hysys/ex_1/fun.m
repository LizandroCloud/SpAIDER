function [ f, df ] = fun( x0, ns, ts, ws, optODE )
global hy
global hysys
% Calculate function values only
 if nargout == 1
 % Forward state integration
 ks=1;
 z0 = [ x0 ];
  for ks = 1:ns
%  [tspan,zs] = ode15s( @(t,x)state(t,x,ws,ks), [ts(ks),ts(ks+1)], z0, optODE );
% hysys = hyconnect;
h = get(hy.activeDocument.Flowsheet.Operations,'Item','SPRDSHT-1');  %
c = hycell(h, {'B1', 'B2', 'B3', 'B4', 'B5'});
hysolvertoggle(hysys);
%hyissolving(hysys);
hyset(c{3},ws(ks));
% hyset(c{1},0);
hysolvertoggle(hy);
%hyintegtoggle(hysys);
%hyissolving(hysys);
hyisintegrating(hy);
%hyintegtime(hysys);
% hyrelease(hysys);
% z0=hyvalue(c(3:5));
  end;
 % Functions
 f = hyvalue(c(3:5))

 % Calculate both function and gradient values
 else
 % Forward state & sensitivity integration
 z0 = [ x0; zeros(3*ns,1) ];
 for ks = 1:ns
 [tspan,zs] = ode15s( @(t,x)sens(t,x,ws,ks), [ts(ks),ts(ks+1)], z0, optODE );
 z0 = zs(end,:)';
 end
 % Functions & Gradients
 f = zs(end,1:3)';
 df = [];
 for is = 1:ns
 df(1:3,is) = zs(end,is*3+1:is*3+3);
 end
 end
 end