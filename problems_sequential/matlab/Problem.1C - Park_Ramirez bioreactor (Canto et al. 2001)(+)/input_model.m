function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study C#
% Problem:     Park Ramirez Bio Reator 
% References:  Eva Canto (2001)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 15; % Final time for integration
input.x0 = [0 0 1 5 1]; %state Initial Condition
input.u_lb = [0]; % control vector lower boundary 
input.u_ub = [2]; % control vector upper boundary 
input.u0 = [0.8]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 5;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 10;  % Final Iteration 
input.x0_lb = [0 0 0 4 0.8]; 
input.x0_ub = [10 10 10 10 10]; 
%--------------------------------------------------------------------------


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