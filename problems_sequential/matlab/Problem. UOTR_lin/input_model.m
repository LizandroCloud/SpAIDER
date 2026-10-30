function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%
%                              #Case study UOTR#
% Problem:     Optimal control of the UOTR model  
% References:  Trica (2022) 
%
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


global PP
% PP=blackoil;
global modelo_u0 modelo_x0 modelo_par
% load('uotr_ctrl_ss_v1.mat','x0','u0','par');
% modelo_x0=x0;
% modelo_u0=u0;
% modelo_par=par;
% global modelo_u0 modelo_x0 modelo_par
% clear x0 u0 par

% idu=[1,5,7,9,15,16,19,20,25];

input.t0 = 0; % Initial time for integration
input.tf = 3600; % Final time for integration
input.x0 = modelo_x0'; %state Initial Condition
% input.u_lb = 1e6/(modelo_u0(2)*(1-modelo_u0(3)))*1e-3; % control vector lower boundary 
input.u_lb = [0.2*modelo_u0(1)*1e-3 1e3];
input.u_ub = [1.1*modelo_u0(1)*1e-3 3.6e3];  % control vector upper boundary 
input.u0 = [modelo_u0(1)*1e-3 3.6e3]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 2;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF =7;  % Final Iteration 

input.x0_lb = modelo_x0*0; %x0 lb
input.x0_ub = abs(modelo_x0)*10.1;  %x0 lb

% algebraic equations
% y(1) = {'(10*x(2)-x(1))'};
%y(1) = {''};

% ode equations
dx(1) = {'u(1)*(10*x(2)-x(1))'};
dx(2) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(3) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(4) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(5) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(6) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(7) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(8) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(9) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(10) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(11) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(12) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(13) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(14) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(15) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(16) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(17) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(18) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(19) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(20) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(21) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(22) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(23) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(24) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(25) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(26) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(27) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};
dx(28) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};

% objective function
J = {'-1+x(end,1)+x(end,2)'};

% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'};
cr_deq(1) = {'[]'};


end