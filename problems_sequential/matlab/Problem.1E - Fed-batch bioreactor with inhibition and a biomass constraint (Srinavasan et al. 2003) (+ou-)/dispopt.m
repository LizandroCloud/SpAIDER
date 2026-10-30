function [ureg1 , ureg0 ] =  dispopt( x0, ns, ts, ws, optODE,iu )


%  cb_in = 5;   % (mol/l)
%  Vo = 1;
%  cao = 0.72;
%  cbo = 0.05;    
%  cd_max = 0.15;
%  cb_max = 0.025;
% Forward state integration; store optimal state & control
z0 = [ x0 ];
topt = [];
xopt = [];
uopt = [];
uopt_lin = [];
tf = 250;
t0 = 0;
reg=0;
ks=1;

     tspan = linspace(t0,tf,500); % t0 = 0; tf = 10, nsteps = 100; 
     [tspan,zs]  = ode45( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE );
    
    topt = [ topt; tspan ];
    xopt = [ xopt; zs ];
    uopt = [ uopt; ws(ks)*ones(length(tspan),1) ];
  
   % perfil nao regularizado... 
   for i=1:(length(topt))
        uopt(i) = u_reg0(topt(i),ts,ws,ks);
   end 
   
    % impressao de perfil regularizado...  
    for i=1:(length(topt))
        u(i) = u_reg1(topt(i),ts,ws,ks);
        u=[u];
    end

   ureg1 = u;
   ureg0 = uopt;
   
   fi=ts(1:ns);

% Mostrando o custo ótimo
% disp(sprintf('Optimal cost for ns=%d: %d', ns));

% Plotando resultados
% figure(1+530)
% plot(topt,xopt(:,1),'r')
% xlabel('t')
% ylabel('Ca(t)')
% title(['ns=',int2str(ns)])
% saveas(gcf,['x1(t)',int2str(ns),'stg.eps'], 'psc2')
% 
% figure(2+530)
% plot(topt,xopt(:,2),'r')
% xlabel('t')
% ylabel('Cb(t)')
% title(['ns=',int2str(ns)])
% saveas(gcf,['x2(t)',int2str(ns),'stg.eps'], 'psc2')

figure(iu+630);
plot(topt,u,'r')
xlabel('t')
ylabel('F(t)')
title(['ns=',int2str(ns)])
% saveas(gcf,['u(t)_',int2str(ns),'stg.eps'], 'psc2')
hold on;

figure(iu+530);
plot(topt,uopt,'b')
xlabel('t')
ylabel('F(t)')
title(['ns=',int2str(ns)])
%saveas(gcf,['u(t)_',int2str(ns),'stg.eps'], 'psc2')
hold on;

figure(iu+530);
plot(fi,ws,'o');
xlabel('t')
ylabel('u(t)')
title(['ns=',int2str(ns)])
% saveas(gcf,['u(t)_',int2str(ns),'stg.eps'], 'psc2')
hold on;

% figure(iu+100+530);
% plot(topt,cc,'r');
% xlabel('t')
% ylabel('cc')
% title(['ns=',int2str(ns)])
% saveas(gcf,['u(t)_',int2str(ns),'stg.eps'], 'psc2')


end
