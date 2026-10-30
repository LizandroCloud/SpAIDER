function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study 1R#
% Problem - Constrained van der Pol Oscillator 
% References:  Vassiliadis et al. (1994) 
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters
par(1)={''}; % L/mol h


% algebraic variables
a_var(1) = {''}; 


% state variables
e_var(1) = {'x(1)'}; e_var(2) = {'x(2)'}; e_var(3) = {'x(3)'};


% control variables
c_var(1) = {'u(1)'};


% ode equations
dx(1) = {'(1-x(2)^2)*x(1)-x(2)+u(1)'};
dx(2) = {'x(1)'};
dx(3) = {'x(1)^2+x(2)^2+u(1)^2'};


% objective function
J = {'x(end,3)'};   % objective functional

% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'-x(:,1)-0.4'} ; % inequality


end