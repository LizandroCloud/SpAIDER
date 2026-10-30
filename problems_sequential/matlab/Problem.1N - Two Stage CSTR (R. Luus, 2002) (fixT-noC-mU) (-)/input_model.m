function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1.N#
% Problem:     Two Stage CSTR
% References:  R. Luus. Iterative dynamic programming. Chapman and Hall/CRC, 2002.
%              
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 2; % Final time for integration
input.x0 = [ 0.15 -0.03 0.10 0 0]; %state Initial Condition
input.u_lb = [-0.5 0.5]; % control vector lower boundary 
input.u_ub = [-0.5 0.5]; % control vector upper boundary 
input.u0 = [0 0]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 
%--------------------------------------------------------------------------


% parameters

par(1)={'tau=0.1'}; %


% algebraic variables

y(1) = {' 0.5 - x(1) - y(3) '};
y(2) = {' -2*( x(2)+0.25 )-u(1)*( x(2) + 0.25 ) + y(3) '};
y(3) = {' ( x(1)+0.5 )*exp( 25 * x(2)*( x(2)+2 )^(-1) )'};
y(4) = {' ( x(3)+0.25 )*exp( 25 * x(4)*( x(4)+2 )^(-1) )'};


% ode equations
dx(1) = {'y(1)'};
dx(2) = {'y(2)'};
dx(3) = {'x(1) - x(3) - tau * y(1) - y(4) + 0.25'};
dx(4) = {'x(2)-2*x(4)-u(2)*( x(4)+0.25 )-tau*(y(2))+y(4)-0.25'};
dx(5) = {'x(1)^2+x(2)^(2)+x(3)^(2)+x(4)^(2)+0.1*(u(1)^(2)+u(2)^(2))'};

% objective function
J = {'x(end,5)'};

% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality


end