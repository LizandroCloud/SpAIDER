function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1.A#
% Problem:SemiBacth Ishothermal CSTR Reactor  
% References:  Srinivasan et al. (2003a) 
%              Schlegel et al. (2005)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 900; % Final time for integration
input.x0 = [0 0]; %state Initial Condition
input.u_lb = [298]; % control vector lower boundary 
input.u_ub = [500]; % control vector upper boundary 
input.u0 = [350]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 


input.x0_lb = [0  ]; %x0 lb
input.x0_ub = [0.8 ]; %x0 ub

%--------------------------------------------------------------------------

% parameters
par(1)={'R = 8.314'}; %
par(2)={'n = 2'}; %
par(3)={'Ea = 63.717*1e3'}; %
par(4)={'A = 7.5581'}; %

% ode equations
dx(1) = {'(1/0.005)*(exp(-Ea/(R*u(1))*(1-x(1))^(n)))'};
dx(2) = {'(u(1))'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-x(end,1) + 1e-4*x(end,2)'};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
% cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end