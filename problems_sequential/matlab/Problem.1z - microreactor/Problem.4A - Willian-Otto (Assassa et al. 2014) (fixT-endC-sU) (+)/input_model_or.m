function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1.A#
% Problem:  Willian-Otto Reactor  
% References:  Assasa and Marquardt (2014) 
%              
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 1000; % Final time for integration
input.x0 = [ 1 0 0 0 0 0 62 2 0]; %state Initial Condition
input.u_lb = [-0.001 20]; % control vector lower boundary 
input.u_ub = [5.79 102]; % control vector upper boundary 
input.u0 = [5  40]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 7;  % Final Iteration 


input.x0_lb = [0 0 0 ]; %x0 lb
input.x0_ub = [0.8 0.1 1.5]; %x0 ub

%--------------------------------------------------------------------------

% parameters
par(1)={'E1 = -6666.7'}; %
par(2)={'E2 = -8333.3'}; %
par(3)={'E3 = -11111'}; %
par(4)={'k1 = 1.6599e6'}; %
par(5)={'k2 = 7.2117e8'}; %
par(6)={'k3 = 2.6745e12'}; %
par(7)={'dh1 = -263.8'}; %
par(8)={'dh2 = -158.3'}; %
par(9)={'dh3 = -226.3'}; % 4.184
par(10)={'cp = 4.184'}; %
par(11)={'cpx = 5554.1'}; %
par(12)={'cex = 125.91'}; %
par(13)={'Tin = 35'}; %
par(14)={'Aref = 9.2903'}; %
par(15)={'mref = 2.1052'}; %
par(16)={'U = 2.3082e-1'}; %
par(17)={'Ma = 100'}; %
par(18)={'Mb = 100'}; %
par(19)={'Mc = 200'}; %
par(20)={'Me = 200'}; %
par(21)={'Mg = 300'}; %
par(22)={'Mp = 100'}; %

par(23)={'roa = 1'}; %
par(24)={'rob = 0'}; %
par(25)={'roc = 0'}; %
par(26)={'roe = 0'}; %
par(27)={'rog = 0'}; %
par(28)={'rop = 0'}; %
par(29)={'tro = 65'}; %
par(30)={'Vo = 2'}; %
par(31)={'Fscal = 1000'}; %


% algebraic equations
 y(1) = {'exp(E1/(x(7) + 273.15))'};
 y(2) = {'exp(E2/(x(7) + 273.15))'};
 y(3) = {'exp(E3/(x(7) + 273.15))'};

% ode equations
dx(1) = {'(-x(1)*u(1)/(Fscal*x(8)))-k1*y(1)*x(1)*x(2)'};

dx(2) = {'((1-x(2))*u(1)/(Fscal*x(8)))-k1*y(1)*x(1)*x(2)- k2*y(2)*x(2)*x(3)'};

dx(3) = {'(-x(3)*u(1)/(Fscal*x(8))) + (Mc*k1*y(1)*x(1)*x(2)/Ma) - (Mc*k2*y(2)*x(2)*x(3)/Mb) - k3*y(3)*x(3)*x(6)'};

dx(4) = {'(-x(4)*u(1)/(Fscal*x(8))) + (Me*k2*y(2)*x(2)*x(3)/Mb)'};

dx(5) = {'(-x(5)*u(1)/(Fscal*x(8))) + (Mg*k3*y(3)*x(3)*x(6)/Mc)'};

dx(6) = {'(-x(6)*u(1)/(Fscal*x(8))) + k2*y(2)*x(2)*x(3) - (Mp*k3*y(3)*x(3)*x(6)/Mc) '};

dx(7) = {'((u(1)*(Tin - x(7)))/(Fscal*x(8))) - (dh1*k1*y(1)*x(1)*x(2)/cp) - (dh2*k2*y(2)*x(2)*x(3)/cp) - (dh3*k3*y(3)*x(3)*x(6)/cp) - (U*Aref/(cp*mref))*(x(7) - u(2))'};

dx(8) = {'u(1)/Fscal'};

dx(9) = {'-cpx*(k2*y(2)*x(2)*x(3) - 0.5*k3*y(3)*x(3)*x(6))*x(8) - cex*(2*k2*y(2)*x(2)*x(3)*x(8))'};

% dx(4) = {'(u(1)*u(1))'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,9) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
% cr_deq(1) = {'x(end,8) - 5.0'} ; % inequality
cr_deq(3) = {'x(:,7) - 90'} ; % inequality
cr_deq(2) = {'-x(:,7) + 60'} ; % inequality
% cr_deq(1) = {'x(end,8) - 5.0'} ; % inequality
end