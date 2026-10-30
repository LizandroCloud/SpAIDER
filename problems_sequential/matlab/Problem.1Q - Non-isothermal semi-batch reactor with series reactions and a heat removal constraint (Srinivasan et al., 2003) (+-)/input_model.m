function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1Q#
% Problem:Problem.1Q - Non-isothermal semi-batch reactor with series reactions and a heat removal constraint 
% References:  Srinivasan et al. (2003) 
%              QUINTO, T; Secchi A.R.; Bicaia Jr. E.C. (2010)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 0.5; % Final time for integration
input.x0 = [ 10 1.1685 0 1]; %state Initial Condition
input.nU = 2; % Number of control variables...
input.u_lb = [(20+273.15) 0]; % control vector lower boundary 
input.u_ub = [ (50+273.15) 1]; % control vector upper boundary 
input.u0 = [20+273.15 0.9]; % scalar initial estimative of control profile for each variable
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 5;  % Final Iteration 
%--------------------------------------------------------------------------


% parameters
par(1)={'k10  = 4'}; % L/mol h
par(2)={'k20 = 800'}; % L/h
par(3)={'E1 = 6e3'}; % J/mol
par(4)={'R = 8.31'}; % J/mol K
par(5)={'dH1 = -3e4'}; % J/mol
par(6)={'dH2 = -1e4'}; % J/mol
par(7)={'qr_maf=1.5e5'}; % J/h
par(8)={'V_maf=1.1'}; % L
par(9)={'E2 = 2e4'}; % J/mol
par(10)={'cb_in = 20'}; % mol/L



y(1) = {'k10 *exp( -E1/(R*u(1)) )'};  % kinetic parameter
y(2) = {'k20 *exp( -E2/(R*u(1)) )'};  % kinetic parameter


% ode equations
dx(1) = {'-y(1)*x(1)*x(2)-( u(2)/x(4) )*x(1)'};
dx(2) = {'-y(1)*x(1)*x(2) + ( (u(2)/x(4))*(cb_in - x(2)))'};
dx(3) = {'y(1)*x(1)*x(2) - y(2)*x(3) - (u(2)/x(4))*x(3)'};
dx(4) = {'u(2)'};

% objective function
J = {'x(end,3)*x(end,4)'};   % objective functional

% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'-dH1*y(:,1)*x(:,1)*x(:,2)*x(:,4) -dH2*y(:,2)*x(:,3)*x(:,4) - qr_maf'} ; % inequality
cr_deq(2) = {'x(end,4) - V_maf'} ; % inequality

end