function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study B#
% Problem:     Optimal control of inequality constrained DAE systems  
% References:  Bell and Sargent, (2000) 
%
% ************************************************************************

% parameters
par(1)={'1'}; %

% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};


% algebraic variables
a_var(1) = {'y(1)'};

% control variables
c_var(1) = {'u(1)'};

% algebraic equations
y(1) = {'(15*x(2)-x(1))'};
%y(1) = {''};

% ode equations
%dx(1) = {'u(1)*(10*x(2)-x(1))'};
dx(1) = {'u(1)*y(1)'};
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