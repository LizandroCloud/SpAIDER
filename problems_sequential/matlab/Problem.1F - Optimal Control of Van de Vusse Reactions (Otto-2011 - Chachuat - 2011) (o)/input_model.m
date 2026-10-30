function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
%                              #Case study 1F#
% Problem:     Cinética de Van der Vusse 
% References:  Otto (2010) (in Portuguese)
% Description: O objetivo de produção neste estudo é maximizar a concentração do componente
%B no tempo final, que é o produto desejado, reduzindo a produção dos componente nãodesejados
%C e D.
% ************************************************************************

[ dx y par e_var a_var c_var J cr_eq cr_deq] = dae_default;

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 6; % Final time for integration
input.x0 = [5.5 0 0 0 378]; %state Initial Condition
input.u_lb = [25 300]; % control vector lower boundary 
input.u_ub = [50 700]; % control vector upper boundary 
input.u0 = [25 355]; % scalar initial estimative of control profile for each variable
input.nU = 2; % Number of control variables...
input.is0 = 4;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 35;  % Final Iteration 
%--------------------------------------------------------------------------



% parameters

par(1)={'V=3.14'}; %
par(2)={'U=642.18'}; %
par(3)={'A=25.1'}; %
par(4)={'k1=1.287e12'}; %
par(5)={'k2=1.287e12'}; %
par(6)={'k3=9.043e9'}; %
par(7)={'E1=9758.3'}; %
par(8)={'E2=9758.3'}; %
par(9)={'E3=8560'}; %
par(10)={'dH1=4200'}; %
par(11)={'dH2=-11000'}; %
par(12)={'dH3=-41850'}; %
par(13)={'ro=1000'}; %
par(14)={'cp=2.8121'}; %
par(15)={'caf=5.5'}; %
par(16)={'cbf=0'}; %
par(17)={'ccf=0'}; %
par(18)={'cdf=0'}; %
par(19)={'tf=378'}; %

% state variables

e_var(1) = {'x(1)'};
e_var(2) = {'x(2)'};
e_var(3) = {'x(3)'};
e_var(4) = {'x(4)'};
e_var(5) = {'x(5)'};

% algebraic variables
a_var(1) = {'y(1)'}; a_var(2) = {'y(2)'}; a_var(3) = {'y(3)'}; a_var(4) = {'y(4)'}; a_var(5) = {'y(5)'}; a_var(6) = {'y(6)'};

y(1) = {'V/u(1)'};  % tal
y(2) = {'U*A*( x(5) - u(2) ) '};  % q
y(3) = {'k1*exp( -E1/x(5) )'};  % k1
y(4) = {'k2*exp( -E2/x(5) )'};  % k2
y(5) = {'k3*exp( -E3/x(5) )'};  % k3
y(6) = {'-( dH1*y(3) + dH2*y(4) + dH3*y(5) )*V '};  % qr

% control variables
c_var(1) = {'u(1)'}; c_var(2) = {'u(2)'}; 

% model
dx(1) = {'( caf - x(1) / y(1) ) - y(3) * x(1) - y(5) * x(1)^2'};
dx(2) = {'( cbf - x(2) / y(1) ) + y(3) * x(1) - y(4) * x(2)'};
dx(3) = {'( ccf - x(3) / y(1) ) +  y(4) * x(2)'};
dx(4) = {'( cdf - x(4) / y(1) ) +  1/2 * x(3)'};
dx(5) = {'( ( tf - x(5) )/ y(1)) + ( ( y(6) - y(2) )/ro*V*cp ) '};

% objective function
%J = {'-1+x(end,1)+x(end,2)'};

J = {'-x(end,2)'};   % calculo da Fobj no tempo final...

% constraints  

% : for 'path', end for 'final'
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'}; % inequality
end