function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study A#
% Problem:SemiBacth Ishothermal CSTR Reactor  
% References:  Srinivasan et al. (2003a) 
%              Schlegel et al. (2005)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters
par(1)={'k1 = 0.053'}; %
par(2)={'k2 = 0.128'}; %
par(3)={'cb_in = 5'}; %
par(4)={'Vo = 1'}; %
par(5)={'lamb=1^-8'}; %
par(6)={'cao=0.72'}; %
par(7)={'cbo=0.05'}; %
par(8)={'cd_mais = 0.15'}; %
par(9)={'cb_mais = 0.025'}; %

% algebraic variables
a_var(1) = {''};

% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
% e_var(4) = {'x(4)'};

% control variables
c_var(1) = {'u(1)'};

% ode equations
dx(1) = {'-k1*x(1)*x(2)-( u(1)/x(3) )*x(1)'};
dx(2) = {'-k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u(1)/x(3))*(cb_in - x(2)))'};
dx(3) = {'u(1)'};
% dx(4) = {'(u(1) - uref)^2'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-cao*Vo + x(end,1)*x(end,3) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'((0.5/x(end,3))*((x(end,1)+cb_in-x(end,2))*x(end,3) - (cao+cb_in-cbo)*Vo)) - cd_mais'} ; % inequality
cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end