function [ f, df ] = fun2( x0, ns, ts, ws, optODE )
global DEBUG;
DEBUG  = 1; 
global reg;
global est;
k1 = 0.053;  % l/(mol min)
k2 = 0.128; % l/(mol min)
cb_in = 5;   % (mol/l)
Vo = 1;
cao = 0.72;
cbo = 0.05;
reg = 1; % regularizacao...
ks=1;
nsteps = 250;
fix =0;  
est = 1;  % integracao por estagios...
tf = 250;
% dasslc parameters
rtol= 1e-10;       % relative tolerance
atol= 1e-10;       % absolute tolerance
rpar=[];       % optional arguments passed to residual and jacobian functions
%tspan = linspace(0,250,nsteps); % t0 = 0; tf = 10, nsteps = 100;
 % Calculando aspenas a funcao
 if nargout == 1
 % Forward state integration
 z0 = [ x0 ];
 
 %=========================================================================
 if fix == 1   % Utilizando solver de passo fixo....
    if reg==0
        if est==1
               for ks = 1:ns
               [zs] = ode5( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks):0.05:ts(ks+1)], [z0]);
               z0 = zs(end,:)';
               end
        else
    %  for ks = 1:ns
    [zs] = ode5( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks):0.3:ts(ns+1)], [z0]);
    %  z0 = zs(end,:)';
    %  end
        end
    else   % se tem regularizacao....
    %    for ks = 1:ns
    [zs] = ode5( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks):2:ts(ns+1)], [z0]);
    %    z0 = zs(end,:)';
    %  end
    end   
 else   % se nao utiliza passo fixo...(ex ode 45...)
     
    
    if reg==0   % se nao regulariza...
       if est==1 
       for ks = 1:ns
       % [tspan,zs]  = dasslc( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ks+1)], z0, optODE );
       rpar(1:length(ws),1) = ws(:,1);
       rpar(length(ws)+1,1) = ks;
       rpar(length(ws)+2:length(ws)+length(ts)+1,1) = ts';
       rpar(length(ts)+length(ws)+2,1) = reg;
       [zp0,ires] = state1(0,z0,zeros(length(z0),1),rpar'); % somente para EDO's
       [t,zs,zps,outp]=dasslc('state1',[ts(ks):0.0001:ts(ks+1)],z0,zp0',rpar',rtol,atol);
       z0 = zs(end,:)';
       end
       else   % se regulariza...
    %    for ks = 1:ns
       rpar(1:length(ws),1) = ws(:,1);
       rpar(length(ws)+1,1) = ks;
       rpar(length(ws)+2:length(ws)+length(ts)+1,1) = ts';
       rpar(length(ts)+length(ws)+2,1) = reg;
       [zp0,ires] = state1(0,z0,zeros(length(z0),1),rpar'); % somente para EDO's
       [t,zs,zps,outp]=dasslc('state1',[0:0.0001:tf],z0,zp0',rpar',rtol,atol);
          
         %[tspan,zs]  = ode45( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
    %    z0 = zs(end,:)';
    %  end
       end
    else
       rpar(1:length(ws),1) = ws(:,1);
       rpar(length(ws)+1,1) = ks;
       rpar(length(ws)+2:length(ws)+length(ts)+1,1) = ts';
       rpar(length(ts)+length(ws)+2,1) = reg;
       [zp0,ires] = state1(0,z0,zeros(length(z0),1),rpar'); % somente para EDO's
       [t,zs,zps,outp]=dasslc('state1',[0:0.0001:tf],z0,zp0',rpar',rtol,atol);
        %[tspan,zs]  = dasslc( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE );
       end
     end
   y=zs;  
%   if ( DEBUG > 0 )
%   if ( outp >= 0 )
% %     disp(sprintf('Integration completed successfully, %d points',length(t)));
%     
%     % result vector contains t in row 1, y(k) in row k
%     for i = 1:size(y,1)
% %       buf = sprintf('% .2e ',y(i,:));
% %       disp(sprintf('pt %3d: t=% .2e y = %s',i,t(i),buf));
%     end
% %     plot(t,y);
%   else
% %     disp(sprintf('Integration failed: result flag = %d',outp));
%   end
% end       
         
 
 % Functions
 f = zs(end,:)';
 cc = (1/f(3))*(cao*Vo - f(1)*f(3));  % calculo de cc no tempo final...
 f(4) = -cc*f(3);   % calculo da Fobj no tempo final...
 
 %Funcao e gradiente..
 else
 % Equacao de sensibilidade...
  zdf = ones(3*ns,1);
  zdf = 0*zdf;
 z0 = [ x0; zdf ];
%  for ks = 1:ns
%  [tspan,zs] = ode45( @(t,x)sens(t,x,ws,ks), [ts(ks),ts(ks+1)], z0, optODE );
%  z0 = zs(end,:)';
%  end
 Rz0 = [ x0 ];
 if reg==0
    %  for ks = 1:ns
    [tspan,zs] = ode45( @(t,x)sens(t,x,ws,ks,ts,reg), [ts(ks),ts(ns+1)], z0, optODE );
    %  z0 = zs(end,:)';
    %  end
 else
    %    for ks = 1:ns
    [tspan,zs] = ode45( @(t,x)sens(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE );
    %    z0 = zs(end,:)';
    %  end
 end
 % Gradientes tb....
 f = zs(end,1:3)';
 cc = (1/f(3))*(cao*Vo - f(1)*f(3));
 f(4) = -cc*f(3);

 
 for is = 1:ns
 df(1:3,is) = zs(end,is*3+1:is*3+3);
 end
 end
 end