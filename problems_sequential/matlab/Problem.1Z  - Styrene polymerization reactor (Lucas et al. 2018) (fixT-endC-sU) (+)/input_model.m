function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1.A#
% Problem:SemiBacth Ishothermal CSTR Reactor  
% References:  Srinivasan et al. (2003a) 
%              Schlegel et al. (2005)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 25; % Final time for integration
input.x0 = [ 6.6832e-2 3.3245 323.56 305.17 2.7547e-4 16.110 1e4*1.3624 0 0]; %state Initial Condition
input.u_lb = [300 0]; % control vector lower boundary 
input.u_ub = [420 100]; % control vector upper boundary 
input.u0 = [471.6]; % scalar initial estimative of control profile for each variable
input.nU = 1; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 6;  % Final Iteration 


input.x0_lb = [0 0 0 ]; %x0 lb
input.x0_ub = [0.8 0.1 1.5]; %x0 ub

%--------------------------------------------------------------------------

% parameters
par(1)={'Ad = 2.1428*10^17'}; %
par(2)={'Ed = 14897'}; %
par(3)={'Ap = 3.816*10^10'}; %
par(4)={'Ep = 3557'}; %
par(5)={'At = 4.50*10^12'}; %
par(6)={'Et = 843'}; %
par(7)={'fi = 0.6'}; %
par(8)={'DeltaHr = -6.99*10^4'}; %
par(9)={'hA = 1.05*10^6'}; %
par(10)={'rhoCp = 1506'}; %
par(11)={'rhocCpc = 4043'}; %
par(12)={'Mm = 104.14'}; %
par(13)={'Qim = 108'}; %
par(14)={'Qi = 108'}; %
par(15)={'Qmm = 378'}; %
par(16)={'Qc =1'}; %
par(17)={'V = 3000'}; %
par(18)={'Vc = 3312.4'}; %
par(19)={'If = 0.5888'}; %
par(20)={'Mf = 8.6981'}; %
par(21)={'Tf = 330'}; %
par(22)={'Tcf = 295'}; %


% algebraic equations

y(1)={'Ad*exp(-Ed/x(3))'};
y(2)={'Ap*exp(-Ep/x(3))'};
y(3)={'At*exp(-Et/x(3))'};
y(4)={'(Qmm/Qim)*Qi '};
y(5)={'1.5*y(4) - Qi'};
y(6)={'Qi + y(5) + y(4)'};
y(7)={'( ( ( 2*fi*y(1)*x(1) )/y(3) )^(0.5) )'};


% ode equations
dx(1) = {'((Qi*If - y(6).*x(1))/V) - y(1).*x(1)'};
dx(2) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};
dx(3) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};
dx(4) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};
dx(5) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};
dx(6) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};
dx(7) = {'((Qi*If - Qt.*x(1))/V) - Kd.*x(1)'};

dx(1) = {'((Qi*If - y(6)*x(1)) / V) - y(1)*x(1)'};
dx(2) = {'( y(4)*Mf - y(6)*x(2) )/ V - y(2)*x(2)*y(7)'};
dx(3) = {'y(6)*(Tf - x(3))/V + (-DeltaHr/rhoCp)*y(2)*x(2)*y(7) -(hA/(rhoCp*V))*(x(3) - x(4))'};
dx(4) = {'u(1)*(Tcf - x(4))/Vc + (hA/(rhocCpc*Vc))*(x(3) - x(4))'};
dx(5) = {'0.5*y(3)*y(7)^2 - (y(6)*x(5))/V'};
dx(6) = {'Mm*y(2)*x(2)*y(7) - (y(6)*x(6))/V'};
dx(7) = {'5*Mm*y(2)*x(2)*y(7) + 3*Mm*((y(2)^2)/y(3))*(x(2)^2) - (y(6)/V)*x(7) '};
dx(8) = {'x(6)*2.5*Qmm*Qi/Qim'};
dx(9) = {'u(1)-320'};

% 
% dx(2) = {'(Qm*Mf - Qt.*x(2))/V - Kp.*x(2).*P'};
% dx(3) = {'Qt.*(Tf - x(3))/V + (-DeltaHr/rhoCp)*Kp.*x(2)*P -(hA/(rhoCp*V)).*(x(3) - x(4))'};
% dx(4) = {'Qc.*(Tcf - x(4))/Vc + (hA/(rhocCpc*Vc)).*(x(3) - x(4))'};
% dx(5) = {'0.5*Kt.*P^2 - (Qt*x(5))/V'};
% dx(6) = {'Mm*Kp.*x(2).*P - (Qt*x(6))/V'};
% dx(7) = {'5*Mm*Kp*x(2)*P + 3*Mm*(Kp^2/Kt)*(x(2)^2) - (Qt/V)*x(7)'};
% dx(8) = {'x(6)*2.5*Qmm*Qi/Qim'};


% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'- x(end,8) '};   % calculo da Fobj no tempo final...
% constraints
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'x(:,3)-387'} ; % inequality
cr_deq(2) = {'104.14*(x(:,7)*x(:,5)/(x(:,6)^2))-1.51'} ; % inequality

end