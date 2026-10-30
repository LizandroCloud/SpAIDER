function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Problem M
% Linear Bang Bang Control
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

% parameters

par(1)={''}; %

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 1; % Final time for integration
input.x0 = [1 0 0]; %state Initial Condition
input.u_lb = [-10 -10]; % control vector lower boundary 
input.u_ub = [10 10]; % control vector upper boundary 
input.u0 = [0 0]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 
%--------------------------------------------------------------------------


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