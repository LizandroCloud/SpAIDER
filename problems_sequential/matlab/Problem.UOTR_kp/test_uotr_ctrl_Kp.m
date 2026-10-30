%Script to run uotr_ctrl_Kp Level-1 MATLAB S-Function
clear, close all

%Physical properties
global PP
PP=blackoil;

%----- Dynamic simulation:
% load('uotr_ctrl_Kp_ss_v0.mat'); u=u0;
% load('uotr_ctrl_Kp_ss_v1.mat'); u=u0;
load('uotr_ctrl_Kp_ss_v1.mat'); u=u0;

%Steady-state #v0
% u(1)=24e3;  %feed.Qf_in
% u(16)=3;    %LIC_001.KC
% u(17)=3;    %LIC_003.KC
% u(18)=3;    %LIC_005.KC

% %Steady-state #v1
% u(1)=18e3;  %feed.Qf_in
% u(16)=4;    %LIC_001.KC
% u(17)=4;    %LIC_003.KC
% u(18)=4;    %LIC_005.KC

% %Steady-state #v2
u(1)=12e3;  %feed.Qf_in
u(16)=5;    %LIC_001.KC
u(17)=5;    %LIC_003.KC
u(18)=5;    %LIC_005.KC

%Jacobian matrix
Nx=size(x0,1);
fhand=@(t,x)uotr_ctrl(t,x,u0,1,[],par);
fval=fhand(0,x0);
[dFdx,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));

JP=abs(dFdx)>=eps; %JPattern
tstart=0; tend=24*3600; tstep=60; tspan=(tstart:tstep:tend)';
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on','JPattern',JP);

tic;
[t,x]=ode15s(@(t,x)uotr_ctrl_Kp(t,x,u,1,[],par),tspan,x0,opt);
timer=toc

figure(1); plot(t,x(:,end));                %[Vg,flare]
figure(2); plot(t,x(:,[6;10;16;27]));       %[SG_001.P; V_TO_001.P; V_TO_002.P; V_001.P]
figure(3); plot(t,x(:,[7;11;17;28]));       %[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
figure(4); plot(t,x(:,[4;8;14]));           %[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
figure(5); plot(t,x(:,[1;12;18]));          %[SG_001.hWCS; TO_001.hW; TO_002.hW]