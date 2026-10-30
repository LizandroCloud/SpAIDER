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
input.tf = 4; % Final time for integration
input.x0 = [ 0 5 1]; %state Initial Condition
input.u_lb = [-40]; % control vector lower boundary 
input.u_ub = [40]; % control vector upper boundary 
input.u0 = [0]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 
input.x0_lb = [0 -5 -30];
input.x0_ub = [400 20 30];
%--------------------------------------------------------------------------

% parameters


% ode equations
dx(1) = {'x(2)^2+0.001*u(1)^2'};
dx(2) = {'x(3)'};
dx(3) = {'u(1)'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,1)'};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
%cr_eq(1) = {'x(end,1)-0'}; % equality
%cr_eq(2) = {'x(end,2)-0'}; % equality
%cr_eq(3) = {'x(end,3)-0'}; % equality
%cr_eq(4) = {'x(end,4)-0'}; % equality
cr_eq(1) = {'[]'} ; % inequality
cr_deq(1) = {'-(x(:,2)-15+(t_out(time)-4)^4) - 0.0'} ; % inequality

end