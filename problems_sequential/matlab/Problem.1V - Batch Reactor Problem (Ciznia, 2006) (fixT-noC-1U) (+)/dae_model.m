function [ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_model

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


optSPAIDER.prob.t0 = 0; % Initial time for integration
optSPAIDER.prob.tf = 1; % Final time for integration
optSPAIDER.prob.x0 = [ 1 0]; %state Initial Condition
optSPAIDER.prob.u_lb =[298]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [398]; % control vector upper boundary 
optSPAIDER.prob.u0 =  [350]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 1; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF = 36;  % Final Iteration 


% algebraic variables
a_var(1) = {'y(1)'};
a_var(2) = {'y(2)'};

y(1) = {'4000*exp(-2500/u(1))'};
y(2) = {'620000*exp(-5000/u(1))'};

% control variables
c_var(1) = {'u(1)'};

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};


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