function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study J#
% Problem:Pressure-constrained batch reactor  
% References:  HUANG et al. (2002) 
%           
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters
par(1)={'k1 = 0.8'}; %
par(2)={'k2 = 0.02'}; %
par(3)={'k3 = 0.003'}; %
par(4)={'N = 1'}; %
par(5)={'R = 8.314'}; %
par(6)={'T = 400'}; %
par(7)={'V = 1'}; %
par(8)={'Pm = 340'}; %

% algebraic variables
a_var(1) = {'y(1)'}; a_var(2) = {'y(2)'};
y(1) = {'V*( x(1)+x(2)+x(3) )'};
y(2) = {'y(1)*R*T/V'};


% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
% e_var(4) = {'x(4)'};

% control variables
c_var(1) = {'u(1)'};

% ode equations
dx(1) = {'-k1*x(1)+k2*x(2)^2 + ( u(1)/V ) - k3*x(1)*x(2)'};
dx(2) = {'k1*x(1) - k2*x(2)^2  - k3*x(1)*x(2)'};
dx(3) = {'k3*x(1)*x(2)'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,3) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'R*T*( x(:,1)+x(:,2)+x(:,3) ) - Pm'} ; % inequality
%cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end