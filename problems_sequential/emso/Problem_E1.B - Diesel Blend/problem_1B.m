% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: []
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% Professors: Argimiro Secchi and Evaristo Biscaia

% ************************************************************************
%                              #Problem E1.B#
% Case: - Diesel Blend
% References: Dyego dos Santos, 2012. (in portuguese)
%             OTIMIZAÇÃO DINÂMICA DE UM PROCESSO DE MISTURA EM LINHA DE DIESEL
% ************************************************************************
                         
                        % Main SpAIDER Interface


%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-6),...
 'TolFun', 10^(-6), 'TolCon', 10^(-6), 'MaxFunEval', 15000000, 'MaxIter', 5e10,...
 'OutputFcn', @myoutput);
%--------------------------------------------------------------------------


% ode parameters
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10);

%--------------------------------------------------------------------------
% method of solution
in_SPAIDER.problem.method = 'single shooting'; 
                          % 'multiple shooting'  (advisory to be paralell if EMSO)
                          % 'simultaneous'
                          % 'indirect'
in_SPAIDER.problem.type = 'optimal control';    
                        % 'dynamic parameter estimation'
                        % 'multiple dynamic objectives'
                        % 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys';   %(If optimset.GradObj is on and/or optimset.GradConstr is on)
                                        % ADA
in_SPAIDER.problem.algorithm = 'parallel'; %serial
%                              'parallel'    % works only for multiple shooting and external routines (as EMSO and Hysys)    
in_SPAIDER.problem.tfree = 0;  % if is free time or fixed time
in_SPAIDER.problem.interpolation = 'constant';  % interpolation of control profile
                                  % 'spline' 
                                  % 'linear'
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 2e-3; % tolerance for adaptation algorithm;
in_SPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
in_SPAIDER.conv.t_ref = 542; % time for CPU reference;
in_SPAIDER.conv.fref =  2.9560; % reference objective functional;
in_SPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'emso'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.est = 1; % option for integration for stages (elements)
in_SPAIDER.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_SPAIDER.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_SPAIDER.casestudy.status = 'off'; % to computate a case study of various thresholds
in_SPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'off';   % 1 for plot or zero for no plot.
in_SPAIDER.plot.state = 0; % for plot state
in_SPAIDER.plot.plot_after = 'off'; % plotting storage .mat
in_SPAIDER.plot.tf_p = 5; % Final time for integration for plotting
in_SPAIDER.plot.record = 'off'; % Record the generated figures  
in_SPAIDER.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'sqtwolog'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes')
in_SPAIDER.wavelets.sorh  = 's';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % on - adaptation wavelets; off -> equidistant refinement...
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.08; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_SPAIDER.wavelets.refinement_position = 0.44; % [0.1 0.9] location of refinement in each stage 
%--------------------------------------------------------------------------

[out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
