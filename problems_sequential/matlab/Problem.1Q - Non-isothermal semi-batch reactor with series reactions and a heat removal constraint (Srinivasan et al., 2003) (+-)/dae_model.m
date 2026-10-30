function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study 1Q#
% Problem:Problem.1Q - Non-isothermal semi-batch reactor with series reactions and a heat removal constraint 
% References:  Srinivasan et al. (2003) 
%              QUINTO, T; Secchi A.R.; Bicaia Jr. E.C. (2010)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


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

% algebraic variables
a_var(1) = {'y(1)'}; a_var(2) = {'y(2)'};

y(1) = {'k10 *exp( -E1/(R*u(1)) )'};  % kinetic parameter
y(2) = {'k20 *exp( -E2/(R*u(1)) )'};  % kinetic parameter

% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};

% control variables
c_var(1) = {'u(1)'};
c_var(2) = {'u(2)'};

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