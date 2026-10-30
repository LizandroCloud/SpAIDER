function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 2V#
% Problem: Instable Problem: Singular Control with Direct Transcription
% References:  Biegler (2010) Book
% *************************************************************************
% 

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = pi/2; % Final time for integration
input.x0 = [ 0 1 0]; %state Initial Condition
input.u_lb = [-1]; % control vector lower boundary 
input.u_ub = [1]; % control vector upper boundary 
input.u0 = [0]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 2;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 

input.x0_lb = [-10 -10 -10]; %x0 lb
input.x0_ub = [10 10 10]; %x0 lb

%--------------------------------------------------------------------------

% parameters


% ode equations
dx(1) = {'x(2)'};
dx(2) = {'u(1)'};
dx(3) = {'0.5*x(2)^2-0.5*x(1)^2'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,3)'};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
%cr_eq(1) = {'x(end,1)-0'}; % equality
%cr_eq(2) = {'x(end,2)-0'}; % equality
%cr_eq(3) = {'x(end,3)-0'}; % equality
%cr_eq(4) = {'x(end,4)-0'}; % equality
cr_eq(1) = {'[]'} ; % inequality
cr_deq(1) = {'[]'} ; % inequality

end