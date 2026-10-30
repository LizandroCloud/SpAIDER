%Script to run dynamic optimization for uotr_ctrl Level-1 MATLAB S-Function
clear, close all

%Physical properties
global PP
PP=blackoil;

%global t ut x Nd Ndx modelo_x0 modelo_u0 ufinal modelo_par modelo_dFdx modelo_dFdu
%global x

%Load initial conditions and model parameters
load('uotr_ctrl_ss.mat','x0','u0','par');

%Jacobians
Nx=size(x0,1);
fhand=@(t,x)uotr_ctrl(t,x,u0,1,[],par); %fhand=@(t,x)uotr_ctrl(t,x,[0;0],[u0(1);u0(1)]*1e-3,u0,1,[],par);
fval=fhand(0,x0);
[modelo_dFdx,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));

Nu=size(u0(1),1);
fhand=@(t,u)uotr_ctrl(t,x0,u,1,[],par); %fhand=@(t,u)uotr_ctrl(t,x0,[0;0],u,u0,1,[],par);
fval=fhand(0,u0); %fval=fhand(0,u0(1)*1e-3);
[modelo_dFdu,~]=numjac(fhand,0,u0,fval,zeros(Nu,1)); %[modelo_dFdu,~]=numjac(fhand,0,u0(1),fval,zeros(Nu,1));

modelo_x0=x0;
modelo_u0=u0;
modelo_par=par;
clear x0 u0 par

%Dynamic optimization parameters
Nd=25;	    %Number of discretization points for uopt
%Ndx=180;     %Number of discretization points for x in csntr evaluation (use only for modelo_eq impulse response)
tstart=0; tstep=60; tend=7200;  t=(tstart:tstep:tend)';	%time span (prediction time)
thrz=1500; %thrz=t(end); %control time horizon
ut=(t(1):(thrz-t(1))/Nd:thrz)'; %Inputs nodal points
%ut=(t(1):(t(end)-t(1))/Nd:t(end))'; 
%xt=(t(1):(t(end)-t(1))/Ndx:t(end))';	%States nodal points (use only for modelo_eq impulse response)
xt=t;
ufinal=1e6/(modelo_u0(2)*(1-modelo_u0(3)));	%ufinal to achieve 1000 Mm^3/d for produced gas

%Guess #1: ramp decrease initial guess with fixed terminal condition for u
%uopt0=(modelo_u0(1):-(modelo_u0(1)-ufinal)/Nd:ufinal)'; uopt0=uopt0(1:end-1);
%ufinal=uopt0(end);

%Guess #2: constant value with variable terminal condition for u
uopt0=modelo_u0(1)*ones(Nd,1);
%ufinal=uopt0(end);

%redim of decision variable for better optimization performance
uopt0=uopt0*1e-3; %redim of decision variable for better optimization performance
ufinal=ufinal*1e-3;

%Plotting initial guess
%ts0=timeseries([uopt0;ufinal],ut); ts0=setinterpmethod(ts0,'zoh');
ts0=timeseries([uopt0;ufinal*ones(1+10,1)],[ut(1:end-1);(thrz:(t(end)-thrz)/10:t(end))']); ts0=setinterpmethod(ts0,'zoh');
figure(1), plot(ts0,'bo-');

%Lower and upper bounds
%uLB=(1-1e-5)*0.1*uopt0(1)*ones(Nd,1); uUB=(1+1e-5)*uopt0(1)*ones(Nd,1);
uLB=(1-1e-5)*ufinal*ones(Nd,1); uUB=(1+1e-5)*uopt0(1)*ones(Nd,1);

%Convoking dynamic optimization
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
tic
%handle_fobj=@(u)objfunc(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu(:,1)); %for modelo_eq impulse response
handle_fobj=@(u)objfunc(u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu(:,1)); %for modelo_eq step response
%handle_fobj=@(u)objfunc(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_par); %for modelo_eq non-linear model
handle_cnstr=@(u)cnstr(t,xt,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu(:,1));
%[cG_fobj,err] = checkGradients(handle_fobj,uopt0,{'Display','on','IsConstraint','false','Tolerance',1e-5});
%[cG_cnstr,err] = checkGradients(handle_cnstr,uopt0,{'Display','on','IsConstraint','true','Tolerance',1e-5});
uopt=fmincon(handle_fobj,uopt0,[],[],[],[],uLB,uUB,handle_cnstr,opt);
timer=toc;

%Plotting optimal input result
%ts=timeseries([uopt;ufinal],ut); ts=setinterpmethod(ts,'zoh');
ts=timeseries([uopt;ufinal*ones(1+10,1)],[ut(1:end-1);(thrz:(t(end)-thrz)/10:t(end))']); ts=setinterpmethod(ts,'zoh');
figure(1); hold on; plot(ts,'ko-'); title('input'); legend('guess','optimal');
%}

%Plotting optimal output result
[xopt,~] = dynmdl(t,uopt,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu(:,1));
%opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on');
%[~,xopt]=ode15s(@(t,x)uotr_ctrl(t,x,ut,[uopt;ufinal]*1e3,modelo_u0,1,[],modelo_par),t,modelo_x0,opt);
figure(2); plot(t,xopt(:,end)); %[SG_001.P; V_TO_001.P; V_TO_002.P; V_001.P]
figure(3); plot(t,xopt(:,[6;10;16;27])); %[SG_001.P; V_TO_001.P; V_TO_002.P; V_001.P]
figure(4); plot(t,xopt(:,[7;11;17;28])); %[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
figure(5); plot(t,xopt(:,[4;8;14])); %[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
figure(6); plot(t,xopt(:,[1;12;18])); %[SG_001.hWCS; TO_001.hW; TO_002.hW]