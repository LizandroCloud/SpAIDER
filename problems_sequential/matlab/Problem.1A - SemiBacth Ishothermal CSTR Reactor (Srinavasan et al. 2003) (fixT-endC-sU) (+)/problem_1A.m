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
% NLP solver options: see help fmincon / optimoptions...
optNLP = optimoptions( 'fmincon', 'Algorithm', 'sqp', ...
 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
 'CheckGradients', false, 'Display', 'iter', 'StepTolerance', 10^(-6),...
 'OptimalityTolerance', 10^(-6), 'ConstraintTolerance', 10^(-6), ...
 'MaxFunctionEvaluations', 1e7, 'MaxIterations', 1e5);

% ODE integration options: see help odeset...
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-5, 'AbsTol', 1e-8);

%--------------------------------------------------------------------------
% method of solution
in_SPAIDER.problem.method = 'single shooting';
                          % 'multiple shooting' (advisable in parallel if EMSO) | 'simultaneous' | 'indirect'
in_SPAIDER.problem.type = 'optimal control';
                          % 'dynamic parameter estimation' | 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys'; % used if SpecifyObjectiveGradient and/or SpecifyConstraintGradient is true
                          % 'psens' | 'automatic diferentiation' | 'numeric'
in_SPAIDER.problem.algorithm = 'serial';
                          % 'parallel' (suggested only for multiple shooting and external routines such as EMSO and Hysys)
in_SPAIDER.problem.tfree = 0;  % 1 for free final time, 0 for fixed final time
in_SPAIDER.problem.interpolation = 'stage';  % interpolation of control profile
                          % 'spline' | 'linear' | 'nearest' | 'stage' | 'reg_const'
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 1e-10; % tolerance of the wavelet adaptation (refinement) loop
in_SPAIDER.conv.f_tol2 = 10;  % second-order tolerance of the adaptation loop
in_SPAIDER.conv.t_ref = 3072.36967484694; % time for CPU reference (s); TODO: record machine/date of this measurement
in_SPAIDER.conv.fref = 0.0; % reference objective functional; TODO: record source (paper or own run)
in_SPAIDER.conv.utarg =0; % 1: use a previous solution as target profile (see spaider_main); 0: default reference
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'ode45'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.reg = 0; % regularization flag; overwritten in spaider_main (reg = 1 only when interpolation is 'reg_const')
in_SPAIDER.solver.sampling = 16; % number of sampling points at each stage
in_SPAIDER.solver.speed_factor = 0; % try to accelerate the solution (may degrade solution quality)
in_SPAIDER.casestudy.status = 'off'; % 'on' runs a case study over several thresholds
in_SPAIDER.casestudy.iter = '64'; % number of case-study runs (string, converted with str2double)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'on';   % 'on' to plot the control profile, 'off' for no plot
in_SPAIDER.plot.state = 0; % 1 to plot the states, 0 for no plot
in_SPAIDER.plot.plot_after = 'off'; % 'on' to plot from stored .mat results instead of running the optimization
in_SPAIDER.plot.tf_p = 250; % final time used for integration when plotting
in_SPAIDER.plot.record = 'off'; % record the generated figures
in_SPAIDER.plot.video = 'off'; % record a video of the control profile refinement
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db2';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'cvpshrink'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_SPAIDER.wavelets.sorh  = 'h';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'mln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % 'on': wavelet-based adaptation; 2: equidistant refinement
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.8; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % maximum number of stages (control discretization)
in_SPAIDER.wavelets.refinement_position = 0.44; % location of refinement in each stage, in [0.1 0.9]
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
