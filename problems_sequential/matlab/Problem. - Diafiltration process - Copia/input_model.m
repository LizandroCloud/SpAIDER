function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study Diafiltration#
%
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 6; % Final time for integration
input.x0 = [ 15 300 0.03]; %state Initial Condition
input.u_lb = [0]; % control vector lower boundary 
input.u_ub = [1]; % control vector upper boundary 
input.u0 = [0.5]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 33;  % Final Iteration 


input.x0_lb = [0 0 0.01 ]; %x0 lb
input.x0_ub = [1600 600 0.035]; %x0 lb

%--------------------------------------------------------------------------

% parameters
par(1)={'s1 = 68.1250*10^-9'}; %
par(2)={'s2 = -56.4512*10^-6'}; %
par(3)={'s3 = 32.5553*10^-3'}; %
par(4)={'s4 = -4.3529*10^-9'}; %
par(5)={'s5 = 3.3216*10^-6'}; %
par(6)={'s6 = -2.7141*10^-3'}; %
par(7)={'z1 = 7.8407*10^-6 '}; %
par(8)={'z2 =  -4.0507*10^-3'}; %
par(9)={'z3 = 1.0585'}; %
par(10)={'z4 =  1.2318*10^-9'}; %
par(11)={'z5 = -9.7660*10^-6'}; %
par(12)={'z6 = -1.1677*10^-3'}; %
par(13)={'w1 = -0.0769*10^-6'}; %
par(14)={'w2 =  -0.0035*10^-3'}; %
par(15)={'w3 = 0.0349*10^-3'}; %
par(16)={'w4 =  0.9961'}; %


% algebraic equations

y(1) = {' ( s1*x(2)^2 + s2*x(2) + s3 )*exp( (s4*x(2)^2 + s5*x(2) + s6) * x(1))  '}; %q
y(2) = {' ( z1*x(2)^2 + z2*x(2) + z3 )*exp( (z4*x(2)^2 + z5*x(2) + z6) * x(1))  '}; %R2
y(3) = {' ( w1*x(2) + w2 )*x(1) + (w3*x(2)+w4) '};  %R1

% ode equations
dx(1) = {'( x(1)/x(3) )*( y(1)*(y(3) - u(1) ) )'};
dx(2) = {'( x(2)/x(3) )*( y(1)*(y(2) - u(1) ) )'};
dx(3) = {' y(1)*(u(1) - 1)'};
% dx(4) = {'(u(1) - uref)^2'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,2)'};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {' x(end,3) - 0.01 '}; % equality
cr_deq(1) = {'-x(:,3) + 0.01 '} ; % inequality
cr_deq(2) = {'x(:,3) - 0.035 '} ; % inequality

end