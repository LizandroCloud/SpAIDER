function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%             #MODELO PSEUDO-HOMOGÊNEO HIDRÓLISE AVICEL#
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 1; % Final time for integration
input.x0 = [ 80 0 1.937 0.1]; %state Initial Condition
input.u_lb = [0 0 100]; % control vector lower boundary [u1 u2] = [Fenz Fsol]
input.u_ub = [0.01 0.04 72*60]; % control vector upper boundary  [u1 u2] = [Fenz Fsol]
input.u0 = [0.001 0.004 72*60]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 5;  % Final Iteration 


input.x0_lb = [0 0 0 0 ]; %x0 lb
input.x0_ub = [80 0 1.937 0.1]; %x0 ub

%--------------------------------------------------------------------------

% parameters
par(1)={'Fout = 0.0'}; %Vazão de saida do reator (L/min)
par(2)={'Eenz = 1.937'}; % [] enzima entra mixer (g/L)
par(3)={'Csol = 80;'}; % [] celulose entra mixer (g/L)
par(4)={'Gsol = 0'}; % [] glicose entra mixer (g/L)
par(5)={'k = 0.0257458'}; % (min^-1)
par(6)={'Km = 11.6035;'}; %(g/L)
par(7)={'Kig = 5.38839'};  %(g/L)

% algebraic equations

y(1) = {' u(2) + u(1)  '}; % Fin
y(2) = {' (u(2)/y(1))*Gsol  '}; %Gin
y(3) = {' (u(1)/y(1))*Eenz '};  %Ein
y(4) = {' k*x(3)*x(1)/(Km*(1+ x(2)/Kig) + x(1))'};  %r %(g/L*min) 
y(5) = {' (u(2)/y(1))*Csol '}; % Cin

% ode equations
dx(1) = {'(y(1)/x(4))*(y(5)-x(1)) - y(4)'}; % C
dx(2) = {'(y(1)/x(4))*(y(2)-x(2)) + 1.111*y(4)'}; % G
dx(3) = {'(y(1)/x(4))*(y(3)-x(3))'}; %  E
dx(4) = {'(y(1)-Fout)'};
% dx(4) = {'(y(1)-Fout)'}; % V



% objective function

 J = {'-x(end,2)*x(end,4)'};
%J = {'-(100*x(end,2))*x(end,4)/(1.111*(80*0.1 + y(5)*y(1).*t_ult))'};
% constraints

cr_eq(1) = {'[]'}; % equality =
cr_deq(1) = {'x(end,4) - 2'} ; % inequality %<=0


end