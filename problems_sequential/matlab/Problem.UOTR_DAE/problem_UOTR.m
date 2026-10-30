% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrosousa@id.uff.br
% UFF - Universidade Federal Fluminense
% Petrobras

% ************************************************************************
%
%                              #Case study UOTR with DAE#
% Problem:     Optimal control of the UOTR model  
% References:  Trica (2022) 
%
%
% ************************************************************************

global PP
PP=blackoil;

global modelo_u0 modelo_par modelo_x0 modelo_dx0
load('uotr_ctrl_dae_ss.mat','x0','u0','par', 'ydot0');
modelo_x0=x0;modelo_x0(end)=eps;
modelo_u0=u0;
modelo_par=par;
modelo_dx0 = ydot0;
clear x0 u0 par ydot0



%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'sqp', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-4),...
 'TolFun', 10^(-6), 'TolCon', 10^(-2), 'MaxFunEval', 15000000, 'MaxIter', 3e10, 'FinDiffType', 'central');
% Ode options: see help ode45...

%Jacobian matrix
Nx=size(modelo_x0,1);
fhand=@(t,x)uotr_ctrl_dae(t,x,modelo_u0,1,[],modelo_par);
fval=fhand(0,modelo_x0);
[dFdx,~]=numjac(fhand,0,modelo_x0,fval,zeros(Nx,1));

%JPattern
JP=abs(dFdx)>=eps;

%Mass matrix

M=[eye(37) 		zeros(37,5)	zeros(37,1)
   zeros(5,37)	zeros(5)	zeros(5,1)
   zeros(1,37)	zeros(1,5)	1];

% ode parameters
%dynamic solver options
%--------------------------------------------------------------------------
optODE=odeset( ...
	'RelTol',1e-5, ...
	'AbsTol',1e-8, ...
	'Mass',M, ...
	'MassSingular','yes', ...
	'MStateDependence','none', ...
	'InitialSlope',modelo_dx0, ...
	'JPattern',JP);
%--------------------------------------------------------------------------
% optODE = odeset( 'RelTol', 1e-5, 'AbsTol', 1e-6, );
% %--------------------------------------------------------------------------
% method of solution
% in_SPAIDER.problem.method = 'single shooting'; 
in_SPAIDER.problem.method = 'single shooting';  %(advisory to be paralell if EMSO)
                          % 'simultaneous'
                          % 'indirect'
in_SPAIDER.problem.type = 'optimal control';    
                        % 'dynamic parameter estimation'
                        % 'multiple dynamic objectives'
                        % 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys';   %(If optimset.GradObj is on and/or optimset.GradConstr is on)
% in_SPAIDER.problem.analitic_derivarive = 'automatic diferentiation';
                                        % 'numeric';
in_SPAIDER.problem.algorithm = 'serial';
%                              'parallel'    % works only for multiple shooting and external routines (as EMSO and Hysys)    
in_SPAIDER.problem.tfree = 0;  % if is free time or fixed time
   in_SPAIDER.problem.interpolation =  'linear';  % interpolation of control profile
%  in_SPAIDER.problem.interpolation =  'spline' ;
 % in_SPAIDER.problem.interpolation =  'linear';
% in_SPAIDER.problem.interpolation =  'nearest';
% in_SPAIDER.problem.interpolation =  'stage';
%  in_SPAIDER.problem.interpolation =  'reg_const';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 1e-12; % tolerance for adaptation algorithm;
in_SPAIDER.conv.f_tol2 = 1e-10;  % second order tolerance for adaptation algorithm;
in_SPAIDER.conv.t_ref = 499; % time for CPU reference;
in_SPAIDER.conv.fref = 0.048055; % reference objective functional;
in_SPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'ode15s'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.est = 1; % option for integration for stages (elements)
in_SPAIDER.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_SPAIDER.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_SPAIDER.casestudy.status = 'off'; % to computate a case study of various thresholds
in_SPAIDER.casestudy.iter = '16';
in_SPAIDER.solver.sampling = 25;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'on';   % 1 for plot or zero for no plot.
in_SPAIDER.plot.state = 0; % for plot state
in_SPAIDER.plot.plot_after = 'off'; % plotting storage .mat
in_SPAIDER.plot.tf_p = 3600*2; % Final time for integration for plotting
in_SPAIDER.plot.record = 'off'; % Record the generated figures  
in_SPAIDER.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'cvpshrink'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_SPAIDER.wavelets.sorh  = 's';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % on - adaptation wavelets; 2 -> equidistant refinement...
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.08; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_SPAIDER.wavelets.refinement_position = 0.42; % [0.1 0.9] location of refinement in each stage 
in_SPAIDER.wavelets.delt = 1.0e-5*0.1;
%--------------------------------------------------------------------------

[out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
figure(2)
bar(out_SPAIDER.time_ac(1:out_SPAIDER.iter+2))