function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 3.Z#
% Problem:Gas Lifting DO  
% References:  
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 20000; % Final time for integration [s]
input.x0 = [ 4350.1 10951 86038 0 0]; %state Initial Condition
input.u_lb = [1 ]; % control vector lower boundary 
input.u_ub = [7 ]; % control vector upper boundary 
input.u0 = [4]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 8;  % Final Iteration 


input.x0_lb = [0; 0; 0; 0; 0 ]; %x0 lb
input.x0_ub = [2e5 ;2e5; 2e5; 2e5; 2e5]; %x0 ub

%--------------------------------------------------------------------------

% parameters
par(1)={'M 	= 0.028'}; %
par(2)={'Ta = 293'}; %
par(3)={'La	= 230.87'}; %
par(4)={'Va	= 22.3e-3'}; %
par(5)={'ro	= 923.9'}; %
par(6)={'Tw	= 293'}; %
par(7)={'Lw	= 2400'}; %
par(8)={'Lr	= 40'}; %
par(9)={'Ar	= 0.0182'}; %
par(10)={'Aw  = 0.0182'}; %
par(11)={'Civ	= 1.6e-3'}; %
par(12)={'Cpc	= 0.5e-2'}; %
par(13)={'Cr 	= 0.8e-4'}; %
par(14)={'rgo = 0'}; %
par(15)={'wgc = 1.36e4'}; %
par(16)={'Ps	= 6276e3'}; %
par(17)={'Pr	= 18632e3'}; %
par(18)={'vo = 1/ro'}; %
par(19)={'R = 8.314'}; %
par(20)={'g = 10'}; %

% algebraic equations

y(4) = {'(R*Ta/(Va*M)+g*La/Va)*x(1)'}; % Pai
y(5) = {'M*y(4)/(R*Ta)'}; % roai
y(6) = {'(x(2)+x(3))/(Lw*Aw)'}; % rom
y(7) = {'(R*Tw/M)*(x(2)/(Lw*Aw-vo*x(3)))'}; % Pwh
y(8) = {'y(7) + (g/Aw)*(x(2)+x(3))'}; % Pwi
y(9) = {'y(8) + ro*g*Lr'}; % Pwb
y(10) = {'Cr*sqrt(ro*(Pr-y(9)))'}; % wro 
y(11) = {'Cpc*sqrt(y(6)*(y(7)-Ps))*u(1)'};	 % wpc
y(12) = {'Civ*sqrt(y(5)*(y(4)-y(8)))'}; %wiv
y(1) = {'x(2)*y(11)/(x(2)+x(3))'}; % wpg
y(2) = {'x(3)*y(11)/(x(2)+x(3))'}; % wpo
y(3) = {'rgo*y(10)'}; % wrg

% ode equations
dx(1) = {'u(1) - y(12)'};
dx(2) = {'y(12) + y(3) - y(1)'};
dx(3) = {'y(10) - y(2)'};
dx(4) = {'(u(2))'};
dx(5) = {'(u(1))'};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-x(end,3)'};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
cr_deq(1) = {'[]'}; %<=0
cr_eq(1) = {'[]'}; % equality
% cr_deq(1) = {'((0.5/x(end,3))*((x(end,1)+cb_in-x(end,2))*x(end,3) - (cao+cb_in-cbo)*Vo)) - cd_mais'} ; % inequality
% cr_deq(2) = {'x(end,2) - cb_mais'} ; % inequality

end