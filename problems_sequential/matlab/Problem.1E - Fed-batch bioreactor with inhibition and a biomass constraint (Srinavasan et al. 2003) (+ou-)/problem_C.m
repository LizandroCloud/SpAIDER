% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Problem C#
%     Fed-batch bioreactor with inhibition and a biomass constraint  
% References:  Srinivasan et al. (2003a) 
%              Visser, Srinivasan, Palanki, & Bonvin, (2000)
%
% ************************************************************************

% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-8),...
 'TolFun', 10^(-8), 'TolCon', 10^(-8), 'MaxFunEval', 20000000);
% Ode options: see help ode45...1
optODE = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8);
f_tol = 0.001;
% ************************************************************************
%problem parameters
t0 = 0; % Initial time for integration
tf = 9; % Final time for integration
tf_p=8; % Final time for plot
t_ref =4000; % Final time for plot
x0 = [ 1; 0; 0; 2 ]; %state Initial Condition
u_lb = [0]; % control vector lower boundary 
u_ub = [1]; % control vector upper boundary 
u0 = [0.1]; % scalar initial estimative of control profile for each variable
nU = 1; % Number of control variables...
is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
isF = 3;  % Final Iteration 
frozen=0;
solver='ode';  % ,fix_runge, ode (for ode45, ode23s, ode15 etc. ),  dassl, hysys or emso.
plot=1;  % 1 for plot or zero for no plot.
plot_state=1; % for plot state
% ************************************************************************

%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa


wname = 'db1';  % wavelets type
tptr = 'sqtwolog'; % thresholding procedure
sorh = 'h';  % hard or soft threshold option
scal = 'mln'; % scale procedure
opC = 1;  % 1 - Wav; 2 -> ADE...
adap=0; % adaption on(thresholding wavelets) or off (fixed)
% [wopt2, erropt2, iout2] = fmincon( @(tf)main(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,optNLP,optODE), tf0, [], [], [],...
% [], 0, 500, [], optNLP);
[wopt, erropt, grad,hessian,wx,tx,nsx,time,time_ac] = spaider(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,frozen,solver,plot,plot_state,tf_p,f_tol,t_ref,optNLP,optODE);
