% A MATLAB example described in detail in the technical report
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% DOTcvp - Dynamic Optimization Toolbox with CVP approach for          %
%  handling continuous and mixed-integer dynamic optimization problems %
% Copyright (C) 2007-2010                                              %
% Tomas Hirmajer et al., thirmajer@gmail.com                           %
%                                                                      %
% The DOTcvp toolbox is completely free of charge under the creative   %
% commons license. The conditions of the license can be found on the   %
% following web page:                                                  %
% http://creativecommons.org/licenses/by-nc-nd/3.0/                    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clear mex; clear all; close all;
% --------------------------------------------------- %
% Initialization:
% --------------------------------------------------- %
data.name                 = 'BifunctionalCatalystBlend';
data.compiler             = 'None'; %['None'|'FORTRAN']

% --------------------------------------------------- %
% Settings for IVP (ODEs, sensitivities):
% --------------------------------------------------- %
data.odes.Def_FORTRAN     = {}; %this option is needed only for FORTRAN parameters definition, e.g. {'double precision k10, k20, ..'}
data.odes.parameters      = {}; %constant parameters before ODE {'T=300','..}
data.odes.Def_MATLAB      = {}; %this option is needed only for MATLAB parameters definition
data.odes.res(1)          = {'-(0.002918487+(-0.008045787)*u(1)+(0.006749947)*u(1)^2+(-0.001416647)*u(1)^3)*y(1)'};
data.odes.res(2)          = {'(0.002918487+(-0.008045787)*u(1)+(0.006749947)*u(1)^2+(-0.001416647)*u(1)^3)*y(1)-((9.509977+(-35.00994)*u(1)+(42.83329)*u(1)^2+(-17.33333)*u(1)^3)+(26.82093+(-95.56079)*u(1)+(113.0398)*u(1)^2+(-44.29997)*u(1)^3))*y(2)+(208.7241+(-719.8052)*u(1)+(827.7466)*u(1)^2+(-316.6655)*u(1)^3)*y(5)'};
data.odes.res(3)          = {'(9.509977+(-35.00994)*u(1)+(42.83329)*u(1)^2+(-17.33333)*u(1)^3)*y(2)'};
data.odes.res(4)          = {'-(0.01921995+(-0.07945320)*u(1)+(0.1105666)*u(1)^2+(-0.05033333)*u(1)^3)*y(4)+(1.350005+(-6.850027)*u(1)+(12.16671)*u(1)^2+(-6.666689)*u(1)^3)*y(5)'};
data.odes.res(5)          = {'(26.82093+(-95.56079)*u(1)+(113.0398)*u(1)^2+(-44.29997)*u(1)^3)*y(2)+(0.01921995+(-0.07945320)*u(1)+(0.1105666)*u(1)^2+(-0.05033333)*u(1)^3)*y(4)-((208.7241+(-719.8052)*u(1)+(827.7466)*u(1)^2+(-316.6655)*u(1)^3)+(1.350005+(-6.850027)*u(1)+(12.16671)*u(1)^2+(-6.666689)*u(1)^3)+(7.339981+(-25.27328)*u(1)+(29.93329)*u(1)^2+(-11.99999)*u(1)^3)+(-0.3950534+(1.679353)*u(1)+(-1.777829)*u(1)^2+(0.4974987)*u(1)^3))*y(5)+(0.1323596+(-0.4696255)*u(1)+(0.5539323)*u(1)^2+(-0.2166664)*u(1)^3)*y(6)+(-0.00002504665+(0.01005854)*u(1)+(-0.01986696)*u(1)^2+(0.009833470)*u(1)^3)*y(7)'};
data.odes.res(6)          = {'(7.339981+(-25.27328)*u(1)+(29.93329)*u(1)^2+(-11.99999)*u(1)^3)*y(5)-(0.1323596+(-0.4696255)*u(1)+(0.5539323)*u(1)^2+(-0.2166664)*u(1)^3)*y(6)'};
data.odes.res(7)          = {'(-0.3950534+(1.679353)*u(1)+(-1.777829)*u(1)^2+(0.4974987)*u(1)^3)*y(5)-(-0.00002504665+(0.01005854)*u(1)+(-0.01986696)*u(1)^2+(0.009833470)*u(1)^3)*y(7)'};
data.odes.black_box       = {'None','1','FunctionName'}; %['None'|'Full'],[penalty coefficient for all constraints],[a black box model function name]
data.odes.ic              = [1.0 0.0 0.0 0.0 0.0 0.0 0.0];
data.odes.NUMs            = size(data.odes.res,2); %number of state variables (y)
data.odes.t0              = 0.0; %initial time
data.odes.tf              = 2000.0; %final time
data.odes.NonlinearSolver = 'Newton'; %['Newton'|'Functional'] /Newton for stiff problems; Functional for non-stiff problems
data.odes.LinearSolver    = 'Dense'; %direct ['Dense'|'Diag'|'Band']; iterative ['GMRES'|'BiCGStab'|'TFQMR'] /for the Newton NLS
data.odes.LMM             = 'Adams'; %['Adams'|'BDF'] /Adams for non-stiff problems; BDF for stiff problems
data.odes.MaxNumStep      = 500; %maximum number of steps
data.odes.RelTol          = 1*10^(-7); %IVP relative tolerance level
data.odes.AbsTol          = 1*10^(-7); %IVP absolute tolerance level
data.sens.SensAbsTol      = 1*10^(-7); %absolute tolerance for sensitivity variables
data.sens.SensMethod      = 'Staggered'; %['Staggered'|'Staggered1'|'Simultaneous']
data.sens.SensErrorControl= 'on'; %['on'|'off']

% --------------------------------------------------- %
% NLP definition:
% --------------------------------------------------- %
data.nlp.RHO              = 10; %number of time intervals
data.nlp.problem          = 'max'; %['min'|'max']
data.nlp.J0               = 'y(7)'; %cost function: min-max(cost function)
data.nlp.u0               = [0.9]; %initial value for control values
data.nlp.lb               = [0.6]; %lower bounds for control values
data.nlp.ub               = [0.9]; %upper bounds for control values
data.nlp.p0               = []; %initial values for time-independent parameters
data.nlp.lbp              = []; %lower bounds for time-independent parameters
data.nlp.ubp              = []; %upper bounds for time-independent parameters
data.nlp.solver           = 'ACOMI'; %['FMINCON'|'IPOPT'|'SRES'|'DE'|'ACOMI'|'MISQP'|'MITS']
data.nlp.SolverSettings   = 'None'; %insert the name of the file that contains settings for NLP solver, if does not exist use ['None']
data.nlp.NLPtol           = 1*10^(-5); %NLP tolerance level
data.nlp.GradMethod       = 'None'; %['SensitivityEq'|'FiniteDifference'|'None']
data.nlp.MaxIter          = 1000; %Maximum number of iterations
data.nlp.MaxCPUTime       = 60*60*0.25; %Maximum CPU time of the optimization (60*60*0.25) = 15 minutes
data.nlp.approximation    = 'PWC'; %['PWC'|'PWL'] PWL only for: FMINCON & without the free time problem
data.nlp.FreeTime         = 'off'; %['on'|'off'] set 'on' if free time is considered
data.nlp.t0Time           = [data.odes.tf/data.nlp.RHO]; %initial size of the time intervals, e.g. [data.odes.tf/data.nlp.RHO] or for the each time interval separately [dt1 dt2 dt3]
data.nlp.lbTime           = 0.01; %lower bound of the time intervals
data.nlp.ubTime           = data.odes.tf; %upper bound of the time intervals
data.nlp.NUMc             = size(data.nlp.u0,2); %number of control variables (u)
data.nlp.NUMi             = 0; %number of integer variables (u) taken from the last control variables, if not equal to 0 you need to use some MINLP solver ['ACOMI'|'MISQP'|'MITS']
data.nlp.NUMp             = size(data.nlp.p0,2); %number of time-independent parameters (p)

% --------------------------------------------------- %
% Equality constraints (ECs):
% --------------------------------------------------- %
data.nlp.eq.status        = 'off'; %['on'|'off'] ECs
data.nlp.eq.NEC           = 1; %number of active ECs
data.nlp.eq.eq(1)         = {''};
data.nlp.eq.time(1)       = data.nlp.RHO;
data.nlp.eq.PenaltyFun    = 'off'; %['on'|'off'] ECs penalty function
data.nlp.eq.PenaltyCoe    = [1.0]; %J0=J0+data.nlp.eq.PenaltyCoe*ViolationOfEqualityConstraint /* only for stochastic solvers */

% --------------------------------------------------- %
% Inequality /path/ constraints (INECs):
% --------------------------------------------------- %
data.nlp.ineq.status      = 'off'; %['on'|'off'] INECs
data.nlp.ineq.NEC         = 2; %number of active INECs
data.nlp.ineq.InNUM       = 1; %how many inequality constraints are '>' else '<'
data.nlp.ineq.eq(1)       = {''};
data.nlp.ineq.eq(2)       = {''};
data.nlp.ineq.Tol         = 0.0005; %tolerance level of violation of INECs
data.nlp.ineq.PenaltyFun  = 'off'; %['on'|'off'] INECs penalty function
data.nlp.ineq.PenaltyCoe  = [1.0 1.0]; %J0=J0+data.nlp.ineq.PenaltyCoe*ViolationOfInequalityConstraint /* for every inequality constraint one parameter */

% --------------------------------------------------- %
% Options for setting of the final output:
% --------------------------------------------------- %
data.options.intermediate = 'off'; %['on'|'off'|'silent'] display of the intermediate results
data.options.display      = 'on'; %['on'|'off'] display of the figures
data.options.title        = 'on'; %['on'|'off'] display of the figure title
data.options.state        = 'on'; %['on'|'off'] display of the state trajectory
data.options.control      = 'on'; %['on'|'off'] display of the control trajectory
data.options.ConvergCurve = 'on'; %['on'|'off'] display of the convergence curve
data.options.Pict_Format  = 'eps'; %['eps'|'wmf'|'both'] save figures as
data.options.report       = 'on'; %['on'|'off'] save data in the dat file
data.options.commands     = {''}; %additional commands, e.g. 'figure(1),.. '
data.options.trajectories = data.odes.NUMs; %how many state trajectories will be displayed
data.options.profiler     = 'off'; %['on'|'off']
data.options.multistart   = 1; %set 1 if the multistart is off, otherwise you have to put here some integer value

data.options.action       = 'single-optimization'; %['single-optimization'|'re-optimization'|'hybrid-strategy'|'simulation']

% --------------------------------------------------- %
% Call of the main function (you do not change this!):
% --------------------------------------------------- %
[data]=dotcvp_main(data);