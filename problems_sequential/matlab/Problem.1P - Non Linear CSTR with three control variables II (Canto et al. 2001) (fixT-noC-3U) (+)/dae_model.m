function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Problem 1G#
%                   Non Linear CSTR
% References:  Canto (2001); 
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


% parameters
par(1)={''}; %

% algebraic variables
a_var(1) = {'y(1)'};

% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};
e_var(5) = {'x(5)'};
e_var(6) = {'x(6)'};
e_var(7) = {'x(7)'};
e_var(8) = {'x(8)'};

% control variables
c_var(1) = {'u(1)'};
c_var(2) = {'u(2)'};
c_var(3) = {'u(3)'};


% algebraic equations
y(1) = {'(u(1) +u(2)+6)'};

% ode equations
dx(1) = {'6 - y(1)*x(1) - 17.6*x(1)*x(2) - 23*x(1)*x(6)*u(3)'};
dx(2) = {'u(1) - y(1)*x(2) - 17.6*x(1)*x(2) - 146*x(2)*x(3)'};
dx(3) = {'u(2) - y(1)*x(3)-73*x(2)*x(3)'};
dx(4) = {'-y(1)*x(4)+35.2*x(1)*x(2)-51.3*x(4)*x(5)'};
dx(5) = {'-y(1)*x(5)+219*x(2)*x(3)-51.3*x(4)*x(5)'};
dx(6) = {'-y(1)*x(6)+102.6*x(4)*x(5)-23*x(1)*x(6)*u(3)'};
dx(7) = {'-y(1)*x(7)+46*x(1)*x(6)*u(3)'};
dx(8) = {'5.8*(y(1)*x(1)-6)-3.7*u(1)-4.1*u(2)+y(1)*(23*x(4)+11*x(5)+28*x(6)+35*x(7))-5*u(3)^2 -0.099'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'(-x(end,8))'};   % calculo da Fobj no tempo final...
% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
end