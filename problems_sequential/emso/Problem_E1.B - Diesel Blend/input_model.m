function [ dx y par e_var a_var c_var J cr_eq cr_deq input] = input_model

% ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www.ceopt.com
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************

% ************************************************************************
%                              #Problem E1.B#
% Case: - Diesel Blend
% References: Dyego dos Santos, 2012. (in portuguese)
%             OTIMIZAÇÃO DINÂMICA DE UM PROCESSO DE MISTURA EM LINHA DE DIESEL
% ************************************************************************

%-------------------------------------------------------------------------
% starting variables
  [ dx y par e_var a_var c_var J cr_eq cr_deq ] = dae_default;
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
input.t0 = 0; % Initial time for integration
input.tf = 150; % Final time for integration
input.x0 = [0 13 13 13 13 13 13 1 6600 1800 42 5 0.15 0.90 5 6 0.3]; %state Initial Condition
input.u_lb = [100 100 100 100 100 100]; % control vector lower boundary 
input.u_ub = [400 400 400 400 400 400 ]; % control vector upper boundary 
input.u0 = [210 100 100 100 178.5 133.77]; % scalar initial estimative of control profile for each variable
input.nU = 6; % Number of control variables...
input.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
input.isF = 3;  % Final Iteration 
input.x0_lb = [ 0 13 13 13 13 13 13 1 6600 1800 42 5 0.15 0.90 5 6 0.3]*0; %x0 lb
input.x0_ub = [ 10 13 13 13 13 13 13 1 6600 1800 42 5 0.15 0.90 5 6 0.3]*10; %x0 lb
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% algebraic equations (for each current)
y(1) = {'Pin + ro*9.8*x(1) - Pdrop'}
y(2) = {'fc*cv*sqrt(Pdrop*df/ro)'};
%--------------------------------------------------------------------------



%-------------------------------------------------------------------------
% ode equations
dx(1) = {'-y(2)'};

%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% objective function
J = {'x(end,3)'}; % objective function
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% constraints  ATTENTION: x(end,1) means a final constraint for variable x(1);
%                         x(:,1)   means a path  constraint for variable x(1);
%                    
%
%cr_eq(1) = {'x(end,1)+x(end,2)'};
%cr_eq(2) = {'x(:,1)^(0.5)+x(:,2)'};
%cr_deq(1) = {'x(:,1)-1'}; %<=0
cr_eq(1) = {'[]'}; % equality
cr_deq(1) = {'[]'} ; % inequality
%-------------------------------------------------------------------------

end