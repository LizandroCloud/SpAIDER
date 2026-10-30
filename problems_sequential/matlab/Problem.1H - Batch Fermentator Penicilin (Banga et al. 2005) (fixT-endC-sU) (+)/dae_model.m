function [ dx y par e_var a_var c_var J cr_eq cr_deq  input] = input_model

% ************************************************************************
%                              #Case study H#
% Problem:     Case study I: optimal control of a fed-batch
%                            fermentor for penicillin production  
% References:  BANGA et al. (2005) 
%
% ************************************************************************
[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;
% parameters

par(1)={''}; %

% algebraic variables
a_var(1) = {''};
% a_var(2) = {'y(2)'};

% y(1) = {'0.11*x(3)/ (0.006*x(1)+x(3))'};
% y(2) = {'0.0055*x(3)/(0.0001 + x(3)*(1 + 10*x(3)))'};

% control variables
c_var(1) = {'u(1)'};

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};

% model
dx(1) = {'(0.11*x(3)/ (0.006*x(1)+x(3)))*x(1)-( u(1)*x(1)/(500*x(4)))'};
dx(2) = {'(0.0055*x(3)/(0.0001 + x(3)*(1 + 10*x(3))))*x(1)-0.01*x(2) - u(1)*x(2)/(500*x(4))'};
dx(3) = {'-((0.11*x(3)/ (0.006*x(1)+x(3)))*x(1)/0.47) - ((0.0055*x(3)/(0.0001 + x(3)*(1 + 10*x(3))))*x(1)/1.2) - (x(1)*0.029*x(3)/(0.0001+x(3))) + (u(1)*((1-x(3)/500)/x(4)))'};
dx(4) = {'u(1)/500'};

% objective function
J = {'( -x(end,2)*x(end,4) )'};   % calculo da Fobj no tempo final...
 
% constraints  
% : for 'path', end for 'final'
 cr_deq(1) = {'x(:,1)-40'};
 cr_deq(2) = {'x(:,3)-25'};
 cr_deq(3) = {'x(:,4)-10'};
 cr_deq(4) = {'x(:,1)+0'};
 cr_deq(5) = {'x(:,3)+0'};
 cr_deq(6) = {'x(:,4)+0'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
%  cr_deq(1) = {'[]'} ; % inequality



end