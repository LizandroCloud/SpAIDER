function [ dx y par e_var a_var c_var J cr_eq cr_deq data] = dae_model

% ************************************************************************
%                              #Case study 2F#
% Problem: Singular Control 2 
% References: CHAPMAN & HALL/CRC Monographs and Surveys in Pure and Applied Mathematics 
%             R. Luus. Iterative dynamic programming. Chapman and Hall/CRC, 2002.
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq ] = dae_default;

data.t0 = 0; % Initial time for integration
data.tf = 5; % Final time for integration
data.x0 = [0 1 0]; %state Initial Condition
data.u_lb = [-1]; % control vector lower boundary 
data.u_ub = [1]; % control vector upper boundary 
data.u0 = [-1]; % scalar initial estimative of control profile for each variable
data.nU = 1; % Number of control variables...
data.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
data.isF = 7;  % Final Iteration 




% parameters
par(1)={''}; %


% algebraic variables
a_var(1) = {''};

% state variables
e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};


% control variables
c_var(1) = {'u(1)'};

% ode equations
dx(1) = {'x(2)'};
dx(2) = {'u(1)'};
dx(3) = {'x(1)^2'};



% objective function
J = {'x(end,3)'}; % objective function


% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality


end