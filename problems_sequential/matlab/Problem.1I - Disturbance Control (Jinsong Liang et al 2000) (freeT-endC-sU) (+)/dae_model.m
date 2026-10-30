function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study C#
% Problem:     Fed-batch bioreactor with inhibition and a biomass constraint  
% References:  Srinivasan et al. (2003a) 
%              Visser, Srinivasan, Palanki, & Bonvin, (2000)
%
% ************************************************************************
[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

input.t0 = 0; % Initial time for integration
input.tf = 1; % Final time for integration
input.x0 = [0.0 0.0]; %state Initial Condition
input.u_lb = [-2 10]; % control vector lower boundary 
input.u_ub = [1 40]; % control vector upper boundary 
input.u0 = [0 15]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 15;  % Final Iteration 


% parameters

par(1)={''}; %

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};



% control variables
c_var(1) = {'u(1)'};


% model
dx(1) = {'x(2)'};
dx(2) = {'u(1)'};


% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,1)'};   % calculo da Fobj no tempo final...

% constraints  
% : for 'path', end for 'final'
 cr_eq(1) = {'x(end,1)-300'};
 cr_eq(2) = {'x(end,2)-0'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
% cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'}; % inequality

end