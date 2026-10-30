function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1R#
% Problem - Constrained van der Pol Oscillator 
% References:  Vassiliadis et al. (1994) 
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters
par(1)={''}; % L/mol h


%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 5; % Final time for integration
input.x0 = [ 0 1 0]; %state Initial Condition
input.nU = 1; % Number of control variables...
input.u_lb = [-0.3 ]; % control vector lower boundary 
input.u_ub = [1 ]; % control vector upper boundary 
input.u0 = [0.5 ]; % scalar initial estimative of control profile for each variable
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 
%--------------------------------------------------------------------------


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