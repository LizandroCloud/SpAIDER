function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study B#
% Problem:     Optimal control of inequality constrained DAE systems  
% References:  Bell and Sargent, (2000) 
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

input.t0 = 0; % Initial time for integration
input.tf = 1; % Final time for integration
input.x0 = [1 0]; %state Initial Condition
input.u_lb = [0]; % control vector lower boundary 
input.u_ub = [1]; % control vector upper boundary 
input.u0 = [0.5]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF =15;  % Final Iteration 

input.x0_lb = [0 0]; %x0 lb
input.x0_ub = [1.5 1.5];  %x0 lb

% algebraic equations
% y(1) = {'(10*x(2)-x(1))'};
%y(1) = {''};

% ode equations
%dx(1) = {'u(1)*(10*x(2)-x(1))'};
dx(1) = {'u(1)*(10*x(2)-x(1))'};
dx(2) = {'u(1)*(x(1)-10*x(2))-(1-u(1))*x(2)'};


% objective function
J = {'-1+x(end,1)+x(end,2)'};

% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'};
cr_deq(1) = {'[]'};
end