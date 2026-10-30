function [ dx y par e_var a_var c_var J cr_eq cr_deq] = input_model

% ************************************************************************
%                              #Problem M
% Linear Bang Bang Control
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

% parameters

par(1)={''}; %

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};

% algebraic variables
a_var(1) = {''};

% control variables
c_var(1) = {'u(1)'};
c_var(2) = {'u(2)'};

% algebraic equations
y(1) = {''};

% ode equations
dx(1) = {'u(2)'};
dx(2) = {'-x(1)+u(1)'};
dx(3) = {'-6*x(1)-12*x(2)+3*u(1)+u(2)'};

% objective function
J = {'(x(end,3))'};   % calculo da Fobj no tempo final...

% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
            

end