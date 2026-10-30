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
fhand=@(t,x)uotr_ctrl(t,x,u0,1,[],par);
fval=fhand(0,x0);
[modelo_dFdx,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));

Nu=size(u0,1);
fhand=@(t,u)uotr_ctrl(t,x0,u,1,[],par);
fval=fhand(0,u0);
[modelo_dFdu,~]=numjac(fhand,0,u0,fval,zeros(Nu,1));

modelo_x0=x0;
modelo_u0=u0;
modelo_par=par;
clear x0 u0 par

%Dynamic optimization parameters
Nd=50;	    %Number of discretization points for uopt
Ndx=100;	%Number of discretization points for x in csntr evaluation
tstart=0; tstep=1; tend=7200;  t=(tstart:tstep:tend)';	%time span
ut=(t(1):(t(end)-t(1))/Nd:t(end))';		%Inputs nodal points
xt=(t(1):(t(end)-t(1))/Ndx:t(end))';	%States nodal points

%Guess #1: ramp decrease initial guess with fixed terminal condition for u
%uopt0=(modelo_u0(1):-(modelo_u0(1)-0.1*modelo_u0(1))/Nd:0.1*modelo_u0(1))';
%ufinal=uopt0(end);

%Guess #2: constant value with variable terminal condition for u
uopt0=modelo_u0(1)*ones(Nd,1);
ufinal=uopt0(end);

%redim of decision variable for better optimization performance
uopt0=uopt0*1e-3; %redim of decision variable for better optimization performance
ufinal=ufinal*1e-3;

%Plotting initial guess
ts0=timeseries([uopt0;ufinal],ut); ts0=setinterpmethod(ts0,'zoh'); %for guess 1
%ts0=timeseries([uopt0;uopt0(end)],ut); ts0=setinterpmethod(ts0,'zoh'); %for guess 2
figure(1), plot(ts0,'bo-');

%Lower and upper bounds
%uLB=(1-1e-5)*0.1*uopt0(1)*ones(Nd,1); uUB=(1+1e-5)*uopt0(1)*ones(Nd,1); %for guess #1
uLB=(1-1e-5)*0.5*uopt0(1)*ones(Nd,1); uUB=(1+1e-5)*uopt0(1)*ones(Nd,1); %for guess #2

%Convoking dynamic optimization
opt=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    ...'EnableFeasibilityMode',true,...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-2, ...
	...'DerivativeCheck','on', ...
    'SpecifyObjectiveGradient',true, ...
	'SpecifyConstraintGradient',true, ...
	'Display','iter');%, ...
    ...'UseParallel', true);
tic
handle_fobj=@(u)fobj(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
handle_cnstr=@(u)cnstr(t,xt,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
uopt=fmincon(handle_fobj,uopt0,[],[],[],[],uLB,uUB,handle_cnstr,opt); %for guess #1
%uopt=fmincon(handle_fobj,uopt0,[],[],[],[],uLB,uUB,handle_cnstr,opt); %for guess #2
timer=toc;

%Plotting optimal input result
%ts=timeseries([uopt;ufinal],ut); ts=setinterpmethod(ts,'zoh'); %for guess #1
ts=timeseries([uopt;uopt(end)],ut); ts=setinterpmethod(ts,'zoh'); %for guess #2
figure(1); hold on; plot(ts,'ko-'); title('input'); legend('guess','optimal');
%}
%Plotting optimal output result
[xopt,~] = modelo_eq(t,uopt,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu); %[x,~] = modelo_eq(uopt0);
%figure(2); plot(t,xopt(:,[6;10;16;27])); %[SG_001.P; V_TO_001.P; V_TO_002.P; V_001.P]
%figure(3); plot(t,xopt(:,[7;11;17;28])); %[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
%figure(4); plot(t,xopt(:,[4;8;14])); %[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
%figure(5); plot(t,xopt(:,[1;12;18])); %[SG_001.hWCS; TO_001.hW; TO_002.hW]