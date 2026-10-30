function [ dx y par e_var a_var c_var J cr_eq cr_deq data] = dae_model

%% ************************************************************************
%                              #Problem J#
%                Lee e Ramirez Bioreactor 
% Paper: Dynamic optimization of chemical and biochemical processes (Canto et al., 2001)
%
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

% parameters

par(1)={''}; %

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};
e_var(5) = {'x(5)'};
e_var(6) = {'x(6)'};
e_var(7) = {'x(7)'};
e_var(8) = {'x(8)'};


% algebraic variables
a_var(1) = {''};

% control variables
c_var(1) = {'u(1)'};
c_var(2) = {'u(2)'};

% algebraic equations
y(1) = {''};

% ode equations
dx(1) = {'u(1)+ u(2)'};
dx(2) = {'((x(3)/(14.35+x(3)+(x(3)^2/111.5))*(x(6)+x(7)*0.22/(0.22+x(5))))*x(2)-(u(1)+u(2))*x(2)/x(1))'};
dx(3) = {'(u(1)/x(1)*100-(u(1)+u(2))*x(3)/x(1)-(x(3)/(14.35+x(3)+(x(3)^2/111.5))*(x(6)+x(7)*0.22/(0.22+x(5))))*x(2)/0.51)'};
dx(4) = {'((0.233*x(3)/(14.35+x(3)+(x(3)^2/111.5))*((0.0005+x(5))/(0.022+x(5))))*x(2)-(u(1)+u(2))*x(4)/x(1))'};
dx(5) = {'(u(2)*4/x(1)-(u(1)+u(2))*x(5)/x(1))'};
dx(6) = {'(-(0.09*x(5)/(0.034+x(5)))*x(6))'};
dx(7) = {'((0.09*x(5)/(0.034+x(5)))*(1-x(7)))'};
dx(8) = {' u(2)'};

% objective function
J = {'(-x(end,1)*x(end,4) + 0*x(end,8))'};   % calculo da Fobj no tempo final...

% constraints  
% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality


end