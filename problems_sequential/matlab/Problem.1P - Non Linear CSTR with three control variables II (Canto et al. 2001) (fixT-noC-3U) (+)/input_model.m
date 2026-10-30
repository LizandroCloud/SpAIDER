function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Problem G#
%                   Non Linear CSTR
% References:  Canto (2001); 
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 0.2; % Final time for integration
input.x0 = [ 0.1883 0.2503 0.0467 0.0899 0.1804 0.1394 0.1046 0.000 ]; %state Initial Condition
input.u_lb = [0 0 0 ]; % control vector lower boundary 
input.u_ub =  [20 6 4 ]; % control vector upper boundary 
input.u0 = [10 3 2 ]; % scalar initial estimative of control profile for each variable
input.nU = 3; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 
%--------------------------------------------------------------------------
input.x0_lb = [0 0 0 0 0 0 0 0]; %x0 lb
input.x0_ub = [10 10 10 10 10 10 10 10]; %x0 ub

% parameters
par(1)={''}; %

% algebraic equations
y(1) = {'(u(1) +u(2)+6)'};

% ode equations
dx(1) = {'6 - y(1)*x(1) - 17.6*x(1)*x(2) - 23*x(1)*x(6)*u(3)'};
dx(2) = {'u(1) - y(1)*x(2) - 17.6*x(1)*x(2) - 146*x(2)*x(3)'};
dx(3) = {'u(2) - y(1)*x(3)-73*x(2)*x(3)'};
dx(4) = {'-y(1)*x(4)+35.2*x(1)*x(2)-51.3*x(4)*x(5)'};
dx(5) = {'-y(1)*x(5)+219*x(2)*x(3)-51.3*x(4)*x(5)'};
dx(6) = {'-y(1)*x(6)+102.6*x(4)*x(5)-23*x(1)*x(6)*u(3)'};
dx(7) = {'-y(1)*x(7)+46*x(1)*x(6)*u(3)'};
dx(8) = {'5.8*(y(1)*x(1)-6)-3.7*u(1)-4.1*u(2)+y(1)*(23*x(4)+11*x(5)+28*x(6)+35*x(7))-5*u(3)^2 -0.099'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'(-x(end,8))'};   % calculo da Fobj no tempo final...
% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
end