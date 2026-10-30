function [ dx y par e_var a_var c_var J cr_eq cr_deq data] = input_model

% ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www.ceopt.com
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************

% ************************************************************************
%                              #Case study 2F#
% Problem: Singular Control 2 
% References: CHAPMAN & HALL/CRC Monographs and Surveys in Pure and Applied Mathematics 
%             R. Luus. Iterative dynamic programming. Chapman and Hall/CRC, 2002.
%
% ************************************************************************

%-------------------------------------------------------------------------
% starting variables
  [ dx y par e_var a_var c_var J cr_eq cr_deq ] = dae_default;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% problem parameters
  data.t0 = 0; % Initial time for integration
  data.tf = 5; % Final time for integration
  data.x0 = [0 1 0]; %state Initial Condition
  data.u_lb = [-1]; % control vector lower boundary 
  data.u_ub = [1]; % control vector upper boundary 
  data.u0 = 0; % initial estimative of control profile for each stage
  data.nU = 1; % Number of control variables...
  data.is0 = 2;  % Initial Iteration (Number of stages is ns^is0 )...
  data.isF = 4;  % Final Iteration 
  data.x0_lb = [-1 -1 -1]; %x0 lb
  data.x0_ub = [1 1 1]; %x0 lb
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% ode equations
dx(1) = {'x(2)'};
dx(2) = {'u(1)'};
dx(3) = {'x(1)^2'};
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