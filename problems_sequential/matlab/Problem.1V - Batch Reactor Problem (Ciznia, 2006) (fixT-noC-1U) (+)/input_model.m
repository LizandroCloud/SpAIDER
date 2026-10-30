function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1V#
% Problem:      Batch Reactor Problem
%                            
% References:  Batch reactor with reactions: A -> B -> C.
%	M. Cizniar, M. Fikar, M. A. Lat, MATLAB Dynamic Optimisation Code DYNOPT. User's Guide, Technical
%	Report, KIRP FCHPT STU Bratislava, Slovak Republic, 2006. 
%
% ************************************************************************
[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

% parameters
par(1)={''}; %


input.t0 = 0; % Initial time for integration
input.tf = 1; % Final time for integration
input.x0 = [ 1 0]; %state Initial Condition
input.u_lb =[298]; % control vector lower boundary 
input.u_ub = [398]; % control vector upper boundary 
input.u0 =  [390]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 2;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 


% algebraic variables
y(1) = {'4000*exp(-2500/u(1))'};
y(2) = {'620000*exp(-5000/u(1))'};


% model
dx(1) = {'-y(1)*x(1)^2'};
dx(2) = {'y(1)*x(1)^2-y(2)*x(2)'};


% objective function
J = {'( -x(end,2) )'};   % calculo da Fobj no tempo final...
 
% constraints  
% Use (:) for 'path constraint', and (end) for 'final constraint'

%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
cr_deq(1) = {'[]'}; %<=0
cr_eq(1) = {'[]'}; % equality
%  cr_deq(1) = {'[]'} ; % inequality



end