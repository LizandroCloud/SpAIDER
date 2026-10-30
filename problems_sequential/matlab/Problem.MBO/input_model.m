function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1.A#
% Problem: Black Oil - EMSO  
% References:  Trica D.J. (2003a) 
%              
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 10800; % Final time for integration
input.x0 = ones(1,2); %state Initial Condition
input.u_lb = [0.4]; % control vector lower boundary 
input.u_ub = [0.5]; % control vector upper boundary 
input.u0 = [0.5]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 


input.x0_lb = [0 0 0 ]; %x0 lb
input.x0_ub = [0.8 0.1 1.5]; %x0 ub

%--------------------------------------------------------------------------

% parameters

% ode equations
dx(1) = {'-x(1)+u(1)'};
dx(2) = {'u(1)'};
% dx(4) = {'(u(1)*u(1))'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-x(end,1) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
% cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end