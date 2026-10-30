function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study 1.N#
% Problem:     Two Stage CSTR
% References:  R. Luus. Iterative dynamic programming. Chapman and Hall/CRC, 2002.
%              
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters

par(1)={'tau=0.1'}; %


% algebraic variables
a_var(1) = {'y(1)'}; a_var(2) = {'y(2)'}; a_var(3) = {'y(3)'}; a_var(4) = {'y(4)'};

y(1) = {' 0.5 - x(1) - y(3) '};
y(2) = {' -2*( x(2)+0.25 )-u(1)*( x(2) + 0.25 ) + y(3) '};
y(3) = {' ( x(1)+0.5 )*exp( 25 * x(2)*( x(2)+2 )^(-1) )'};
y(4) = {' ( x(3)+0.5 )*exp( 25 * x(4)*( x(4)+2 )^(-1) )'};


% state variables
e_var(1) = {'x(1)'}; e_var(2) = {'x(2)'};  e_var(3) = {'x(3)'}; e_var(4) = {'x(4)'}; e_var(5) = {'x(5)'};

% control variables
c_var(1) = {'u(1)'}; c_var(2) = {'u(2)'};

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