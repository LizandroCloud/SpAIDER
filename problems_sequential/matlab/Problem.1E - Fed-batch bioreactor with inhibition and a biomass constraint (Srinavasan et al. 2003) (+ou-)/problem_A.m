% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Problem A#
%                    SemiBacth Ishothermal CSTR Reactor
% References:  Schlegel (2004), PhD Thesis; 
% Srinavan et al. (2003);
%
% ************************************************************************

% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'on', 'GradConstr', 'on',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 15000, 'MaxIter', 10);
% Ode options: see help ode45...
optODE = odeset( 'RelTol', 1e-7, 'AbsTol', 1e-7);

% ************************************************************************
%problem parameters
t0 = 0; % Initial time for integration
tf = 250; % Final time for integration
x0 = [ 0.72; 0.05; 1 ]; %state Initial Condition
u_lb = [0]; % control vector lower boundary 
u_ub = [0.001]; % control vector upper boundary 
u0 = [0.0005]; % scalar initial estimative of control profile for each variable
nU = 1; % Number of control variables...
is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
isF = 10;  % Final Iteration 
frozen=0;
solver='ode';  % ,fix_runge, ode (for ode45, ode23s, ode15 etc. ),  dassl, hysys or emso.
plot=1;  % 1 for plot or zero for no plot.
% ************************************************************************

%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa


wname = 'db1';  % wavelets type
tptr = 'rigrsure'; % thresholding procedure
sorh = 'h';  % hard or soft threshold option
scal = 'mln'; % scale procedure
opC = 1;  % 1 - Wav; 2 -> ADE...
adap=0; % adaption on(thresholding wavelets) or off (fixed)
% [wopt2, erropt2, iout2] = fmincon( @(tf)main(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,optNLP,optODE), tf0, [], [], [],...
% [], 0, 500, [], optNLP);
[wopt, erropt, grad,hessian,wx,tx] = main(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,frozen,solver,plot,optNLP,optODE);
