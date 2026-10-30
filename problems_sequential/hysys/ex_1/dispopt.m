function dispopt( x0, ns, ts, ws, optODE, iu )

% Forward state integration; store optimal state & control
z0 = [ x0 ];
topt = [];
xopt = [];
uopt = [];
for ks = 1:ns
[tspan,zs] = ode15s( @(t,x)state(t,x,ws,ks), [ts(ks),ts(ks+1)], z0, optODE );
 z0 = zs(end,:)';
 topt = [ topt; tspan ];
 xopt = [ xopt; zs ];
 uopt = [ uopt; ws(ks)*ones(length(tspan),1) ];
 end

 % Display Optimal Cost
 disp(sprintf('Optimal cost for ns=%d: %d', ns, zs(end,3)));

 % Plot results
 figure(1)
 plot(topt,xopt(:,1),'r',topt,xopt(:,2),'b')
 xlabel('t')
 ylabel('x1(t), x2(t)')
 title(['ns=',int2str(ns)])
 saveas(gcf,['x1(t)_x2(t)_',int2str(ns),'stg.eps'], 'psc2')

 figure(2)
 plot(xopt(:,1),xopt(:,2),'r')
 xlabel('x1(t)')
 ylabel('x2(t)')
 title(['ns=',int2str(ns)])
 saveas(gcf,['x2(x1)_',int2str(ns),'stg.eps'], 'psc2')

 figure(iu+3)
 plot(topt,uopt,'r')
 xlabel('t')
 ylabel('u(t)')
 title(['ns=',int2str(ns)])
 saveas(gcf,['u(t)_',int2str(ns),'stg.eps'], 'psc2')
 end