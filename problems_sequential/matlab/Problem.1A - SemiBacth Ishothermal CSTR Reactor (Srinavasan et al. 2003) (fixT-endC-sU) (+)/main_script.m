% *************************************************************************
% Script Template for SPAIDER toolbox
%
% Author: Lizandro de Sousa Santos
% Contact: lizandrosousa@id.uff.br
% Institution: Universidade Federal Fluminense (UFF) - Brazil
% Department of Chemical Engineering 
%
% Description:
%
% This script is the main template for solving dynamic optimization and
% optimal control problems using the SPAIDER framework. It provides
% a structured configuration of numerical, solver, and refinement
% options, which can be adapted to different case studies.
%
% In this specific instance, the template is instantiated for a
% semi-batch, isothermal Continuous Stirred Tank Reactor (CSTR),
% where the objective is to compute optimal control profiles subject
% to nonlinear dynamic constraints defined by mass balance equations.
%
% The script is intended to serve as a starting point for new SPAIDER
% applications, allowing users to modify the process model, objective
% function, constraints, and numerical settings while preserving a
% consistent solution workflow.
% ------------------------------------------------------------------------
% Methodology:
%
% - Optimal control formulation with fixed final time
% - Single-shooting discretization of control variables
% - Nonlinear programming solved using SQP (fmincon)
% - Adaptive control refinement using wavelet-based techniques
% - External ODE solver interface (ode45 by default)
% ------------------------------------------------------------------------
% Key Features:
%
% - Adaptive wavelet-based refinement of control stages
% - Analytical gradients for objective and constraints
% - Configurable solver, discretization, and convergence parameters
% - Optional plotting and case study analysis
% ------------------------------------------------------------------------
% References:
%
% Srinivasan, B. et al. (2003a)
% Schlegel, M. et al.   (2005)
%
% Notes:
%
% - The script assumes the SPAIDER toolbox is correctly installed
%   and initialized.
% - Numerical tolerances and wavelet thresholds should be adjusted
%   according to the physical scaling of the problem.
% - For stiff reactor models, implicit ODE solvers (e.g., ode15s)
%   may be more appropriate.
%
% *************************************************************************


% ************************************************************************
% Environment initialization
% ************************************************************************

startup
% Runs the startup script. Used to configure paths, load SPAIDER dependencies,
% initialize global variables, and prepare the execution environment.


% ************************************************************************
% Nonlinear optimization (NLP) options — fmincon settings
% ************************************************************************

optNLP = optimset( ...
    'Algorithm', 'sqp', ...        % Optmization Algorithms ('sqp', 'interior-point', 'active-set')
    'GradObj', 'on', ...           % Objective function gradient: Provided analytically (on) or Numerically (off)
    'GradConstr', 'on', ...        % Constraint gradients: Provided analytically (on) or Numerically (off)
    'DerivativeCheck', 'off', ...  % Numerical gradient checking: 'on' or 'off'
    'Display', 'iter', ...         % Displays optimization progress at each iteration
    'TolX', 10^(-6), ...           % Tolerance on the decision variables
    'TolFun', 10^(-6), ...         % Tolerance on objective function changes
    'TolCon', 10^(-6), ...         % Constraint violation tolerance
    'MaxFunEval', 15000000, ...    % Maximum number of function evaluations
    'MaxIter', 10e10);             % Maximum number of iterations (effectively unlimited)

% ************************************************************************
% ODE solver options
% ************************************************************************

optODE = odeset( ...
    'RelTol', 1e-5, ...            % Relative tolerance of the ODE solver
    'AbsTol', 1e-8);               % Absolute tolerance (more stringent)

% ************************************************************************
% Problem definition — solution method
% ************************************************************************

in_SPAIDER.problem.method = 'single shooting';
% Single shooting method:
% – Only the control variables are optimized
% – State variables are obtained by direct integration of the ODEs
% – Simple and efficient, but sensitive for stiff or highly nonlinear systems

% Alternative methods:
% 'multiple shooting'  → more robust, parallelizable
% 'simultaneous'       → states and controls discretized together
% 'indirect'           → based on optimality conditions (Pontryagin)

in_SPAIDER.problem.type = 'optimal control';
% Problem type:
% – Optimal control (default)
% Alternatives:
% – dynamic parameter estimation
% – static optimization

in_SPAIDER.problem.analitic_derivarive = 'psens';
% Analytical derivative computation method:
% 'psens' → parametric sensitivities (efficient and accurate)
% Alternatives:
% 'automatic diferentiation'
% 'numeric'
% 'sensitivy analisys'

in_SPAIDER.problem.algorithm = 'serial';
% Serial execution mode (single processor)
% 'parallel' → recommended for multiple shooting and external simulators

in_SPAIDER.problem.tfree = 0;
% 0 → fixed final time
% 1 → free final time (optimized)

in_SPAIDER.problem.interpolation = 'stage';
% Control profile interpolation:
% 'stage' → piecewise-constant control
% Alternatives:
% 'spline', 'linear', 'nearest', 'reg_const'

% ************************************************************************
% Convergence criteria for the adaptive algorithm
% ************************************************************************

in_SPAIDER.conv.f_tol = 1e-10;
% Primary tolerance for the adaptive refinement algorithm

in_SPAIDER.conv.f_tol2 = 10;
% Secondary (higher-order) tolerance criterion

in_SPAIDER.conv.t_ref = 3072.36967484694;
% Reference CPU time (used for normalization or comparison)

in_SPAIDER.conv.fref = 0.0;
% Reference value of the objective functional

in_SPAIDER.conv.utarg = 0;
% Target control value (used in specific adaptive criteria)

% ************************************************************************
% Integration and solver options
% ************************************************************************

in_SPAIDER.solver.type = 'external';
% 'external' → uses external solvers (MATLAB ODE solvers, Aspen, Hysys)
% 'internal' → uses SPAIDER internal routines

in_SPAIDER.solver.name = 'ode45';
% Selected solver:
% ode45 → explicit Runge–Kutta method (non-stiff problems)
% Alternatives: ode23s, dasslc, external simulators

in_SPAIDER.solver.reg = 0;
% Regularization option (only active when estimation is disabled)

in_SPAIDER.solver.sampling = 16;
% Number of sampling points per stage

in_SPAIDER.solver.speed_factor = 0;
% Heuristic acceleration option (may reduce solution accuracy)

in_SPAIDER.casestudy.status = 'off';
% 'on' → computes a case study for multiple thresholds

in_SPAIDER.casestudy.iter = '64';
% Number of iterations for the case study

% ************************************************************************
% Plotting options
% ************************************************************************

in_SPAIDER.plot.control = 'on';
% Enables control profile plotting

in_SPAIDER.plot.state = 0;
% Enables state trajectory plotting (0 → disabled)

in_SPAIDER.plot.plot_after = 'off';
% Stores plots in .mat files after execution

in_SPAIDER.plot.tf_p = 250;
% Final time used only for plotting purposes

in_SPAIDER.plot.record = 'off';
% Automatically records generated figures

in_SPAIDER.plot.video = 'off';
% Records a video of control refinement over iterations

% ************************************************************************
% Adaptive wavelet-based refinement options
% ************************************************************************

in_SPAIDER.wavelets.wname = 'db2';
% Wavelet family (Daubechies order 2)

in_SPAIDER.wavelets.tptr = 'cvpshrink';
% Thresholding procedure:
% cvpshrink → cross-validation-based shrinkage

in_SPAIDER.wavelets.sorh = 'h';
% 'h' → hard thresholding
% 's' → soft thresholding

in_SPAIDER.wavelets.scal = 'mln';
% Multi-level normalization strategy

in_SPAIDER.wavelets.procedure = 'on';
% Enables wavelet-based adaptive refinement

in_SPAIDER.wavelets.thresholding = 'automatic';
% Automatic threshold selection

in_SPAIDER.wavelets.norm_threshold = 0.8;
% Norm-based threshold parameter

in_SPAIDER.wavelets.fix_threshold = 1.0e-5;
% Fixed threshold value (if selected)

in_SPAIDER.wavelets.max_refinement = 128;
% Maximum number of refinement stages

in_SPAIDER.wavelets.refinement_position = 0.44;
% Relative position within each stage where refinement occurs

in_SPAIDER.wavelets.delt = 1.0e-55 * 0.1;
% Minimum admissible increment (prevents infinite refinement loops)

in_SPAIDER.wavelets.rel = 'on';
% Enables relative refinement strategy

% ************************************************************************
% Main execution with error handling
% ************************************************************************

try
    % Main SPAIDER execution call
    [out_SPAIDER] = spaider_main(in_SPAIDER, optNLP, optODE, SPAIDER);

catch ME
    % Error handling block

    idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
    % Extracts the last segment of the error identifier

    if strcmp(idSegLast,'')
        % Special case: empty identifier segment
        % [out_SPAIDER] = spaider_main(in_SPAIDER, optNLP, optODE);
    else
        % Re-executes without the SPAIDER argument
        [out_SPAIDER] = spaider_main(in_SPAIDER, optNLP, optODE);
    end
end
