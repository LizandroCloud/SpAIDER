function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1E#
% Problem: Feed-batch Bioreactor with inibition and biomass constraint 
% References: Sirinivasan, 2003
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 8; % Final time for integration
input.x0 = [1 0 0 2]; %state Initial Condition
input.u_lb = [0]; % control vector lower boundary 
input.u_ub =  [1]; % control vector upper boundary 
input.u0 =  [0.5]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 20;  % Final Iteration 
input.x0_lb = [0.7 0 0 1.8]; %x0 lb
input.x0_ub = [1 0.1 0.1 3];  %x0 lb
%--------------------------------------------------------------------------

% parameters
par(1)={'mim = 0.53'}; % l/(h)
par(2)={'km = 1.2'}; % l/(h)
par(3)={'ki = 22'};  % (g/L)
par(4)={'Yx = 0.4'};  % (g/L)
par(5)={'Yp = 1'};  % (g/L)
par(6)={'v = 0.5'};  %L/h;
par(7)={'Sin = 20'};  % g/l
par(8)={'xo = 1'};  %g/L;
par(9)={'So = 0'};  %g/L;
par(10)={'Po = 0'};  %g/L;
par(11)={'Vo = 2;'};  %L
par(12)={'xmax= 3'}; %g/L;
par(13)={'smax= 5'}; %g/L;

% algebraic equations
y(1) = {'mim*x(2) / (km + x(2)+ ( (x(2)^2)/ki ))'};


% ode equations

dx(1) = {'y(1)*x(1) - u(1)*x(1)/x(4)'};
dx(2) = {'- y(1)*x(1)/Yx - v*x(1)/Yp + u(1)*(Sin-x(2))/x(4)'};
dx(3) = {' v*x(1) - u(1)*x(3)/x(4)'};
dx(4) = {'u(1)'};


% objective function
J = {'-x(end,3)*x(end,4)'};   % calculo da Fobj no tempo final...

% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_deq(1) = {'x(:,1) - xmax'}; % equality
cr_deq(2) = {'x(:,2) - smax'}; % equality
cr_eq(1) = {'[]'} ; % inequality


end