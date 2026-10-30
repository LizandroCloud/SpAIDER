% ************************************************************************
% Author:      Lizandro de Sousa Santos
% Email:       lizandrosousa@id.uff.br
% Institution: Universidade Federal Fluminense
%
% ************************************************************************
%                              #Case study 1.A#
% Problem: Semi-Batch Isothermal CSTR Reactor
% References:  Srinivasan et al. (2003a)
%              Schlegel et al. (2005)
%
% ************************************************************************
if exist('startup', 'file') ~= 2
    error('problem_1A:startup', 'startup.m not found on the MATLAB path; run from the SpAIDER root.');
end
startup
%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimoptions( 'fmincon', 'Algorithm', 'sqp', ...
 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
 'CheckGradients', false, 'Display', 'iter', 'StepTolerance', 10^(-10),...
 'OptimalityTolerance', 10^(-10), 'ConstraintTolerance', 10^(-3), ...
 'MaxFunctionEvaluations', 1e7, 'MaxIterations', 1e5);
% ConstraintTolerance is looser than the other tolerances on purpose (assumed:
% fref was obtained with it); tighten it only together with a new fref.

% ode parameters: see help odeset...
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-9, 'AbsTol', 1e-9);

%--------------------------------------------------------------------------
% method of solution
in_SPAIDER.problem.method = 'single shooting';
                          % 'multiple shooting' (advisable in parallel if EMSO) | 'simultaneous' | 'indirect'
in_SPAIDER.problem.type = 'optimal control';
                          % 'dynamic parameter estimation' | 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys'; % used if SpecifyObjectiveGradient and/or SpecifyConstraintGradient is true
                          % 'automatic diferentiation' | 'numeric'
in_SPAIDER.problem.algorithm = 'serial';
                          % 'parallel' (suggested only for multiple shooting and external routines such as EMSO and Hysys)
in_SPAIDER.problem.tfree = 0;  % 1 for free final time, 0 for fixed final time
in_SPAIDER.problem.interpolation = 'stage';  % interpolation of control profile
                          % 'spline' | 'linear' | 'nearest' | 'stage' | 'reg_const'
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 1e-5; % tolerance for adaptation algorithm;
in_SPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
in_SPAIDER.conv.t_ref = 3072.36967484694; % time for CPU reference (s); TODO: record machine/date of this measurement
in_SPAIDER.conv.fref = 0.431721; % reference objective functional; TODO: record source (paper or own run)
in_SPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'ode45'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_SPAIDER.solver.sampling = 50; % number of sampling points at each stage
in_SPAIDER.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_SPAIDER.casestudy.status = 'off'; % to computate a case study of various thresholds
in_SPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'on';   % 1 for plot or zero for no plot.
in_SPAIDER.plot.state = 0; % for plot state
in_SPAIDER.plot.plot_after = 'off'; % plotting storage .mat
in_SPAIDER.plot.tf_p = 250; % Final time for integration for plotting
in_SPAIDER.plot.record = 'off'; % Record the generated figures  
in_SPAIDER.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_SPAIDER.wavelets.sorh  = 'h';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % on - adaptation wavelets; 2 -> equidistant refinement...
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.8; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_SPAIDER.wavelets.refinement_position = 0.44; % [0.1 0.9] location of refinement in each stage 
in_SPAIDER.wavelets.delt = 1.0e-56; % offset added to the threshold in spaider_main (negligible)
in_SPAIDER.wavelets.rel = 'on';
%--------------------------------------------------------------------------

try
    [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE,SPAIDER);
catch ME
    % fallback only when the optional SPAIDER argument is missing
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') && contains(ME.message, 'SPAIDER')
        [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE);
    else
        rethrow(ME);
    end
end
