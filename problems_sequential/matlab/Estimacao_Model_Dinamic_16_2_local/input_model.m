function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Mestrado Mayara#
%                       Estimacao_Model_Dinamic_16_2_local
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 0.1; % Final time for integration
input.x0 = [ 1.5 5.3640E-02 3.1337E+01 0 0 0 0 0 0.0072 0.2642]; %state Initial Condition
input.u_lb = [1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4]; % control vector lower boundary 
input.u_ub = [1 0.1 10 10 10 10 10 10 10 10 10]; % control vector upper boundary 
input.u0 = [1 1 1 1 1 1 1 1 1 1 1 ]; % scalar initial estimative of control profile for each variable
input.nU = 11; % Number of control variables...
input.is0 = 1;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 


input.x0_lb = [0 0 0 ]; %x0 lb
input.x0_ub = [0.74 0.07 1.22]; %x0 lb

%--------------------------------------------------------------------------

% parameters

% ode equations
dx(1) = {'-k1*x(1)*x(2)-( u(1)/x(3) )*x(1)'};
dx(2) = {'-k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u(1)/x(3))*(cb_in - x(2)))'};
dx(3) = {'u(1)'};
dx(4) = {'u(1)'};
dx(5) = {'u(1)'};
dx(6) = {'u(1)'};
dx(7) = {'u(1)'};
dx(8) = {'u(1)'};
dx(9) = {'u(1)'};
dx(10) = {'u(1)'};


% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,10) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality


end