% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Case study B#
% Problem:     Optimal control of inequality constrained DAE systems  
% References:  Bell and Sargent, (2000) 
%
% ************************************************************************

% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-11),...
 'TolFun', 10^(-11), 'TolCon', 10^(-11), 'MaxFunEval', 15000000, 'MaxIter', 5e10);

%--------------------------------------------------------------------------
% method of solution
in_SPAIDER.problem.method = 'single shooting'; 
                          % 'multiple shooting'  (advisory to be paralell if EMSO)
                          % 'simultaneous'
                          % 'indirect'
in_SPAIDER.problem.type = 'optimal control';    
                        % 'dynamic parameter estimation'
                       
                        % 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys';   %(If optimset.GradObj is on and/or optimset.GradConstr is on)
% in_SPAIDER.problem.analitic_derivarive = 'automatic diferentiation';
                                        % 'numeric';
in_SPAIDER.problem.algorithm = 'serial';
%                              'parallel'    % suggered for only for multiple shooting and external routines (as EMSO and Hysys)    
in_SPAIDER.problem.tfree = 0;  % if is free time or fixed time
in_SPAIDER.problem.interpolation = 'constant';  % interpolation of control profile
                                  % 'spline' 
                                  % 'linear'
%--------------------------------------------------------------------------


% convergence parameters
% ************************************************************************
optODE = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10);
%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 5e-7; % tolerance for adaptation algorithm;
in_SPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
in_SPAIDER.conv.t_ref = 3072.36967484694; % time for CPU reference;
in_SPAIDER.conv.fref = 0.428601; % reference objective functional;
in_SPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------
% optSPAIDER.conv.t_ref = 2.017877925843673e+03; % time for CPU reference;
% optSPAIDER.conv.fref = 14.693; % reference objective functional;
%--------------------------------------------------------------------------
optSPAIDER.prob.t0 = 0; % Initial time for integration
optSPAIDER.prob.tf = 1; % Final time for integration
optSPAIDER.prob.x0 = [1 0]; %state Initial Condition
optSPAIDER.prob.u_lb = [0]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [1]; % control vector upper boundary 
optSPAIDER.prob.u0 = [0.5]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 1; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF =37;  % Final Iteration 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'ode45'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.est = 0; % option for integration for stages (elements)
in_SPAIDER.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_SPAIDER.solver.sampling = 9; % number of sampling points at each stage
in_SPAIDER.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_SPAIDER.casestudy.status = 'off'; % to computate a case study of various thresholds
in_SPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'on';   % 1 for plot or zero for no plot.
in_SPAIDER.plot.state = 0; % for plot state
in_SPAIDER.plot.plot_after = 'off'; % plotting storage .mat
in_SPAIDER.plot.tf_p = 1; % Final time for integration for plotting
in_SPAIDER.plot.record = 'off'; % Record the generated figures  
in_SPAIDER.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------


% wavelet options
%--------------------------------------------------------------------------
%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa
%  'gcv'for generalized cross validation
%  'bayes'for bayes theorem

% wname = 'db1';  % wavelets type
% tptr = 'heursure'; % thresholding procedure
% sorh = 'h';  % hard or soft threshold option
% scal = 'sln'; % scale procedure
% Wav = 1;  % 1 - Wav; 2 -> ADE...
% adap=2; % adaption on(thresholding wavelets) or off (fixed)

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_SPAIDER.wavelets.sorh  = 's';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % on - adaptation wavelets; 2 -> equidistant refinement...
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.08; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_SPAIDER.wavelets.refinement_position = 0.74; % [0.1 0.9] location of refinement in each stage 
in_SPAIDER.wavelets.delt = 1.0e-5*0.1;
%--------------------------------------------------------------------------

try
[out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE,SPAIDER)
catch ME
    
        idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match'); 
        
         if strcmp(idSegLast,'') 
                 [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
         else
                 [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
         end
   
end
