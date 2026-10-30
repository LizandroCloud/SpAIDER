function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Problem Diafiltration#
% Problem:  
% References:  
%              
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;


input.t0 = 0; % Initial time for integration
input.tf = 8000; % Final time for integration
input.x0 = [ 0 11.227 325 0 ]; %state Initial Condition
input.u_lb = [270]; % control vector lower boundary 
input.u_ub = [400]; % control vector upper boundary 
input.u0 = [300]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 


% parameters
par(1)={'R = 8.314'}; %
par(2)={'MWa = 104'}; %
par(3)={'mwa = 205.7'}; %
par(4)={'Tref = 355.382'}; %
par(5)={'Tamb = 280.382'}; %
par(6)={'V = 0.00567'}; %
par(7)={'B1 = 0.193'}; %
par(8)={'P = 1.594'}; %
par(9)={'B2 = 0.167'}; %
par(10)={'Cppa = 3.140'}; %
par(11)={'Cpma = 1.675'}; %
par(12)={'Cpwa = 4.18'}; %
par(13)={'tb = 4000'}; %
par(14)={'a1_ina = 0.000*7.5*0.001'}; %
% par(13)={'a1_ina = 7.5*0.001'}; %
par(15)={'koa= 55'}; %
par(16)={'k1a = 1000'}; %
par(17)={'k2a = 0.4'}; %
par(18)={'Ea = 29560.888'}; %
par(19)={'coa = 5.2e-5'}; %
par(20)={'c1a = 16.44'}; %
par(21)={'c2a = 2.3'}; %
par(22)={'c3a = 1.563'}; %
par(23)={'aoa = 555.556'}; %
par(24)={'Hpa = -70152.16'}; %
par(25)={'gma = 900'}; %
par(26)={'gpa = 1040'}; %
par(27)={'gwa = 1000'}; %
par(28)={'doa = .814'}; %
par(29)={'d1a = -5.13'}; %
par(30)={'hf = 1/0.7'}; %
   


% algebraic variables

y(1) = {'x(2)/(x(1)+x(2)+ mwa)'};
y(2) = {'(coa*exp(c1a*y(1)))*10^(c2a*(aoa/x(3))-c3a)'};
y(3) = {'koa*(exp(-Ea/(R*x(3))))*(k1a*y(2))^k2a'};
y(4) = {'x(1)/MWa'};
y(5) = {'1.2*y(3)*y(4)'};
y(6) = {'((x(1)/gma) + (x(2)/gpa) + (mwa/gwa))*P/B1 + B2'};
y(7) = {'0.5*(x(3)+u(1))'};
y(8) = {'(coa*exp(c1a*y(1)))*(10^(c2a*(aoa/y(7))-c3a))'};
y(9) = {'doa*exp(d1a*y(8))'};
y(10)= {'1.0/(1.0/y(9) + 1.0/hf)'};
y(11)= {'(( a1_ina*Cpma*(Tamb-x(3)) ) - ( y(9)*y(6)*(x(3)-u(1)) ) - ( V*(x(3)-Tamb) ) - y(5)*Hpa )'};
y(12)= {'( x(1)*Cpma + x(2)*Cppa + mwa*Cpwa )'};

% ode equations
dx(1) = {'a1_ina-y(5)*MWa'};
dx(2) = {'y(5)*MWa'};
dx(3) = {'y(11)/y(12)'};
dx(4) = {'( x(3)-Tref )^2'};
% dx(4) = {'(u(1) - uref)^2'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'x(end,4)'};   % calculo da Fobj no tempo final...
% constraints
cr_eq(1) = {'[]'};
cr_deq(1) = {'[]'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
% cr_eq(1) = {'[]'}; % equality
% cr_deq(1) = {'x(:,3)-0.035'} ; % inequality
% cr_deq(2) = {'x(:,3)+0.01'} ; % inequality
% cr_deq(3) = {'x(end,3)+0.01'} ; % inequality

end