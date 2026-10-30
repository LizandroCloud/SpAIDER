function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

% ************************************************************************
%                              #Case study C#
% Problem:     Park Ramirez Bio Reator 
% References:  Eva Canto (2001)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

% parameters

par(1)={''}; %

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};
e_var(5) = {'x(5)'};

% algebraic variables
a_var(1) = {''};

% control variables
c_var(1) = {'u(1)'};

% model
dx(1) = {'(4.75*(21.88*x(4)/( (x(4) + 0.4)*(x(4)+62.5) ) ) / (0.12 + (21.88*x(4)/( (x(4) + 0.4)*(x(4)+62.5) ) )))*(x(2)-x(1))-u(1)*x(1)/x(5)'};
dx(2) = {'( ( x(4)/(0.1+x(4)) )*exp(-5.01*x(4)) )*x(3) - u(1)*x(2)/x(5)'};
dx(3) = {'(21.88*x(4)/( (x(4) + 0.4)*(x(4)+62.5) ) )*x(3) - u(1)*x(3)/x(5)'};
dx(4) = {'-7.3*(21.88*x(4)/( (x(4) + 0.4)*(x(4)+62.5) ) )*x(3) + u(1)*(20-x(4))/x(5)'};
dx(5) = {'u(1)'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-x(end,1)*x(end,5)'};   % calculo da Fobj no tempo final...

% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'}; % inequality
end