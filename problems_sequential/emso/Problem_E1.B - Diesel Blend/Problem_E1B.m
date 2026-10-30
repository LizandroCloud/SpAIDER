% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: 
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil

% ************************************************************************
%                              #Problem E1.B#
% Case: - Diesel Blend
% References: Dyego dos Santos, 2012. (in portuguese)
%             OTIMIZAÇÃO DINÂMICA DE UM PROCESSO DE MISTURA EM LINHA DE DIESEL
% ************************************************************************

%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-4),...
 'TolFun', 10^(-4), 'TolCon', 10^(-4), 'MaxFunEval', 15000000, 'MaxIter', 5e10);
% Ode options: see help ode45...

% ODE parameters
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-4, 'AbsTol', 1e-4);

% problem parameters
%--------------------------------------------------------------------------
% method of solution
optSPAIDER.problem.type = 'multiple shooting'; 
                       % 'multiple shooting'  (advisory to be paralell if emso)
                       % 'simultaneous'
optSPAIDER.problem.tfree = 0;
% convergence parameters
%--------------------------------------------------------------------------
optSPAIDER.conv.f_tol = 0.00047; % tolerance for adaptation algorithm;
optSPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
optSPAIDER.conv.t_ref = 256.32; % time for CPU reference;
optSPAIDER.conv.fref = 0.2685; % reference objective functional;
optSPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
optSPAIDER.prob.t0 = 0; % Initial time for integration
optSPAIDER.prob.tf = 150*60; % Final time for integration
optSPAIDER.prob.x0 = [0 13 13 13 13 13 13 1 6600 1800 42 5 0.15 0.90 5 6 0.3]; %state Initial Condition
optSPAIDER.prob.u_lb = [100 100 100 100 100 100]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [400 400 400 400 400 400 ]; % control vector upper boundary 
optSPAIDER.prob.u0 = [210 100 100 100 178.5 133.77]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 6; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF = 3;  % Final Iteration 
%--------------------------------------------------------------------------

%integration options

%--------------------------------------------------------------------------
optSPAIDER.solver.name = 'emso';
optSPAIDER.solver.wfact = 0.4;
optSPAIDER.solver.type = 'external';
optSPAIDER.solver.reg = 0;
optSPAIDER.solver.est = 0;
optSPAIDER.solver.frozen = 0;
optSPAIDER.casestudy.state = 0;
optSPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------

% plot options
%--------------------------------------------------------------------------
optSPAIDER.plot.control = 1;   % 1 for plot or zero for no plot.
optSPAIDER.plot.state = 0; % for plot state
optSPAIDER.plot.plot_after = 0; % plotting storage .mat
optSPAIDER.plot.tf_p = 150; % Final time for integration for plotting
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
optSPAIDER.wavelets.tunn = 64; % constraints for number of stages
optSPAIDER.wavelets.wfact = 0.4; % [0.1 0.9] location of refinement in each stage 
%--------------------------------------------------------------------------

[SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)
