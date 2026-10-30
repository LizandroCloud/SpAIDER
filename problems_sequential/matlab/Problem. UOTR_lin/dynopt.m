%Run the lines below in MATLAB(C) prompt to call dynopt for
%uotr_ctrl.m
%{
clear
global PP; PP=blackoil;
global modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu
Nd=25; Ndx=30; tspan=(0:60:7200)'; thrz=1500;
[uopt,timer] = dynopt(Nd,Ndx,tspan,thrz);
%}

function [uopt,timer] = dynopt(Nd,Ndx,tspan,thrz)
%Nd: 	Number of discretization points for uopt
%Ndx: 	Number of discretization points for states in constraint function
%tspan:	Time span nodal points for modelo_eq.m integration
%thrz:	Horizon control time

close all
global modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu modelo_par

%Load initial conditions and model parameters
load('uotr_ctrl_ss.mat','x0','u0','par');

%Jacobians
Nx=size(x0,1);
fhand=@(t,x)uotr_ctrl(t,x,u0,1,[],par); %fhand=@(t,x)uotr_ctrl(t,x,[0;0],[u0(1);u0(1)]*1e-3,u0,1,[],par);
fval=fhand(0,x0);
[modelo_dFdx,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));

Nu=size(u0(1),1);
fhand=@(t,u)uotr_ctrl(t,x0,[u;u0(2:end)],1,[],par); %fhand=@(t,u)uotr_ctrl(t,x0,[0;0],u,u0,1,[],par);
fval=fhand(0,u0(1)); %fval=fhand(0,u0(1)*1e-3);
[modelo_dFdu,~]=numjac(fhand,0,u0(1),fval,zeros(Nu,1)); %[modelo_dFdu,~]=numjac(fhand,0,u0(1),fval,zeros(Nu,1));

modelo_x0=x0;
modelo_u0=u0;
modelo_par=par;
clear x0 u0 par

%Dynamic optimization parameters
ut=(tspan(1):(thrz-tspan(1))/Nd:thrz)'; 				%Inputs nodal points
%xt=(tspan(1):(tspan(end)-tspan(1))/Ndx:tspan(end))'; 	%States nodal points
ufinal=1e6/(modelo_u0(2)*(1-modelo_u0(3)));	%ufinal to achieve 1000 Mm^3/d for produced gas

%Initial guess
uopt0=modelo_u0(1)*ones(Nd+1,1);

%redim of decision variable for better optimization performance
ufinal=ufinal*1e-3; 
uopt0=uopt0*1e-3;

%Plotting initial guess
plot_uopt0=timeseries([uopt0(1:end-1);uopt0(end)*ones(1+10,1)],[ut(1:end-1);(thrz:(tspan(end)-thrz)/10:tspan(end))']);
plot_uopt0=setinterpmethod(plot_uopt0,'zoh');
figure(1), plot(plot_uopt0,'bo-');

%Lower and upper bounds
uLB=(1-1e-5)*ufinal*ones(Nd+1,1); uUB=(1+1e-5)*uopt0(1)*ones(Nd+1,1);

%call fmincon
opt=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'TolX', 1e-6, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-2, ...
	...'DerivativeCheck','on', ...
    'GradObj','on', ...
	'GradConstr','on', ...
	'Display','iter');%, ...
    ...'UseParallel', true);

uLast = []; % Last place computeall was called
myf = [];
myJ = []; % Use for objective at xLast
mygradJ = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint
mygradc = []; % Use for nonlinear inequality constraint
mygradceq = []; % Use for nonlinear equality constraint

%Internal parameters for SPAIDER which must be set
ks=[];
reg=[];
ns=[];
ne=[];
te=[];
n_s=Nx;				%Number of states
n_c=Nu;				%Number of inputs
ns_a=[Nd+1];		%Number of decision points for each input
ts_frozen=[];
points=[];
uws=[]; %uws=ones(Nu,Nd+1);	%Decision variables incidence matrix for sensitivity matrix
tal=1;
est=[];
yi=[];
utarg=[];
tstarg=[];
nstarg=[];

obj_hndl 	= @(u)fun(tspan,u,ks,ut,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
cnstr_hndl	= @(u)cfun(tspan,u,ks,ut,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);

% Call fmincon
tic
uopt=fmincon(obj_hndl,uopt0,[],[],[],[],uLB,uUB,cnstr_hndl,opt);
    function [J,gradJ] = fun(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
        if ~isequal(ws,uLast) % Check if computation is necessary
            [myf,myJ,mygradJ,myc,myceq,mygradc,mygradceq] = ...
                computeall(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
            uLast = ws;
        end
        % Now compute objective function
        f = myf;
        J = myJ;
        gradJ = mygradJ;
    end

    function [c,ceq,gradc,gradceq] = cfun(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
        if ~isequal(ws,uLast) % Check if computation is necessary
            [myf,myJ,mygradJ,myc,myceq,mygradc,mygradceq] = ...
                computeall(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
            uLast = ws;
        end
        % Now compute constraint function
        f = myf;
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
        gradc = mygradc;
        gradceq = mygradceq;
    end
timer=toc;

%Plotting optimal input result
plot_uopt=timeseries([uopt(1:end-1);uopt(end)*ones(1+10,1)],[ut(1:end-1);(thrz:(tspan(end)-thrz)/10:tspan(end))']);
plot_uopt=setinterpmethod(plot_uopt,'zoh');
figure(1); hold on; plot(plot_uopt,'ko-'); title('input'); legend('guess','optimal');

% %Plotting 
% funh_model=@(t,x)modelo_eq(t,x,[uopt;uopt(end)],ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
% [~,f_out]=ode15s(funh_model,tspan,[modelo_x0,odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on'));
% xopt=f(:,1:size(modelo_x0,1));
% figure(2); plot(tspan,xopt(:,end)); 			%[Vg,flare]
% figure(3); plot(tspan,xopt(:,[6;10;16;27])); 	%[SG_001.P; V_TO_001.P; V_TO_002.P; V_001.P]
% figure(4); plot(tspan,xopt(:,[7;11;17;28])); 	%[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
% figure(5); plot(tspan,xopt(:,[4;8;14])); 		%[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
% figure(6); plot(tspan,xopt(:,[1;12;18])); 		%[SG_001.hWCS; TO_001.hW; TO_002.hW]

end