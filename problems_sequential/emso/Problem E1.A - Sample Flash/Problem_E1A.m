% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peresetallq.coppe.ufrj.br
% home page: www.ceopt.com
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil

% ************************************************************************
%                              #Problem E1.A#
% Case: - Sample Flash 
% References: 
%
% ************************************************************************
% matlabpool close;
% poolSize = 4;
% matlabpool (poolSize);
%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 15000000, 'MaxIter', 5e10);

% Ode options: see help ode45...
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-6);

% problem parameters
%--------------------------------------------------------------------------
% method of solution
optSPAIDER.problem.type = 'multiple shooting'; 
                       % 'multiple shooting'  (advisory to be paralell if emso)
                       % 'simultaneous'

% convergence parameters
%--------------------------------------------------------------------------
optSPAIDER.conv.f_tol = 0.00047; % tolerance for adaptation algorithm;
optSPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
optSPAIDER.conv.t_ref = 256.32; % time for CPU reference;
optSPAIDER.conv.fref = 0.2685; % reference objective functional;
optSPAIDER.wavelets.tunn = 64; % constraints for number of stages
optSPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
optSPAIDER.prob.t0 = 0; % Initial time for integration
optSPAIDER.prob.tf = 800; % Final time for integration
optSPAIDER.prob.x0 = [338 0.4 0.1 0.1 0.1 0.1 0.1 0.1]; %state Initial Condition
optSPAIDER.prob.u_lb = [300 0.19*[338 0.4 0.1 0.1 0.1 0.1 0.1 0.1] ]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [700 1.91*[338 0.4 0.1 0.1 0.1 0.1 0.1 0.1] ]; % control vector upper boundary 
optSPAIDER.prob.u0 = [500 338 0.4 0.5 0.1 0.1 0.1 0.1 0.1]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 1; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF = 5;  % Final Iteration 
%--------------------------------------------------------------------------

%integration options

%--------------------------------------------------------------------------
optSPAIDER.solver.name = 'emso';
optSPAIDER.wavelets.wfact = 0.4;
optSPAIDER.solver.type = 'external';
optSPAIDER.solver.reg = 0;
optSPAIDER.solver.est = 1;
optSPAIDER.solver.frozen = 0;
optSPAIDER.problem.tfree = 0;
optSPAIDER.casestudy.state = 0;
optSPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------

% plot options
%--------------------------------------------------------------------------
optSPAIDER.plot.control = 1;   % 1 for plot or zero for no plot.
optSPAIDER.plot.state = 1; % for plot state
optSPAIDER.plot.plot_after = 0; % plotting storage .mat
optSPAIDER.plot.tf_p = 800; % Final time for integration for plotting
%--------------------------------------------------------------------------

% wavelet options
%--------------------------------------------------------------------------
%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa


% wname = 'db1';  % wavelets type
% tptr = 'rigrsure'; % thresholding procedure
% sorh = 'h';  % hard or soft threshold option
% scal = 'sln'; % scale procedure
% Wav = 1;  % 1 - Wav; 2 -> ADE...
% adap=2; % adaption on(thresholding wavelets) or off (fixed)

%--------------------------------------------------------------------------
optSPAIDER.wavelets.wname = 'db1';  % wavelets type
optSPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure
optSPAIDER.wavelets.sorh  = 'h';  % hard or soft threshold option
optSPAIDER.wavelets.scal  = 'sln'; % scale procedure
optSPAIDER.wavelets.Wav   = 2;  % 1 - Wav; 2 -> ADE...
optSPAIDER.wavelets.adap  = 1; % adaption on(thresholding wavelets) or off (fixed)
optSPAIDER.wavelets.normfrac  = 0.08; % norm
optSPAIDER.wavelets.thr  = 1.0e-5; % fixed
%--------------------------------------------------------------------------

[SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)
