% teste integrador dasslc, ode45, fixo (por estagios e direto) 
global DEBUG;
DEBUG  = 1; 
global reg;
global est;
optODE = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10);

nsteps = 250;
fix =1;  
% dasslc parameters
rtol= 1e-10;       % relative tolerance
atol= 1e-10;       % absolute tolerance

rpar=[];       % optional arguments passed to residual and jacobian functions
%tspan = linspace(0,250,nsteps); % t0 = 0; tf = 10, nsteps = 100;

 % Calculando aspenas a funcao
x0 = [ 0.72; 0.05; 1];
t0 = 0.0; % Tempo Inicial
tf = 250; % Tempo Final
ns = 8;

% estagios de controle
%ts = [ 0  9.6154 19.2308 28.8462 38.4615 48.0769 57.6923 67.3077 76.9231 86.5385 96.1538  105.7692 115.3846 125.0001 134.6154 144.2308 153.8462 158.6538 163.4615 168.2692 173.0769 182.6923 192.3077 197.1154 201.9231 206.7308 211.5385 216.3462 221.1538 225.9615 230.7692 240.3846 250.0000];
ts = [0.000000001
31.2500000000001
62.5000000000001
93.7500000000001
125.0001
156.250000000001
187.500000000001
218.750000000001
250.0000001];

% perfil de controle otimo
%ws = 1.0e-003 *[1.0000 1.0000 0.8760 0.8301 0.8165 0.7956 0.7777 0.7601 0.7435 0.7276 0.7125 0.6980 0.6842 0.6710 0.6585 0.6461 0.6365 0.6325 0.6262 0.6224 0.6136 0.6037 0.5890 0.6173 0.4904 0.0000 0.0000 -0.0000 -0.0000 0 -0.0000 0];
ws =[0.00100000000000000
0.000790924655873986
0.000757394745083947
0.000706906978075635
0.000641631177033681
0.000708961644478830
0.000228010099942998
6.47809602783118e-26
];
z0 = [ x0 ];
topt = [];
xopt = [];
uopt = [];
uopt_lin = [];
tf = 250;
t0 = 0;
reg=0;
ks=1;
est = 1;  % integracao por estagios...
tout=[];
zout=[]; 
% tspan = linspace(t0,tf,500); % t0 = 0; tf = 10, nsteps = 100; 
% [tspan,zs]  = ode45( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE );
% topt = [ topt; tspan ];
% xopt = [ xopt; zs ];
% uopt = [ uopt; ws(ks)*ones(length(tspan),1) ];

 if fix == 1   % Utilizando solver de passo fixo....
    if reg==0  % Regularizacao
        if est==1  % Estagios...
               for ks = 1:ns
               tspan = linspace(t0,ts(ks+1),(1000/(ns))); % t0 = 0; tf = 10, nsteps = 100; 
               [zs] = ode5( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks):0.25:ts(ks+1)], [z0]);
               z0 = zs(end,:)';
               tt=ts(ks):0.25:ts(ks+1);
               tout = [ tout; tt' ];
               zout = [ zout; zs ];
               uopt = [ uopt; ws(ks)*ones(length(tspan),1) ];
               end
        else  % Nao tem estagios...
    %  for ks = 1:ns
    tspan = linspace(t0,tf,1001); % t0 = 0; tf = 10, nsteps = 100; 
    [zs] = ode5( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks):0.25:ts(ns+1)], [z0]);
    topt = [ topt; tspan ];
    xopt = [ xopt; zs ];
       for i=1:(length(topt))
        uopt(i) = u_reg0(topt(i),ts,ws,ks);
        end 
%     uopt = [ uopt; ws(ks)*ones(length(tspan),1) ];
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
       tout=[tout; t]; 
       zout=[zout; zs]; 
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
       tout=[tout t]; 
       zout=[zout zs]; 
       end
     end

 figure(1)
 stairs(tout,zout(:,1),'r')
 xlabel('t')
 ylabel('Ca(t)')
 title(['ns=',int2str(ns)])

 
 figure(2)
 plot(tout,zout(:,2),'r')
 xlabel('t')
 ylabel('Cb(t)')
 title(['ns=',int2str(ns)])

