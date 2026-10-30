% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www.ceopt.com
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil

% ************************************************************************
%                              #Case study 2F#
% Problem: Singular Control 2 
% References: CHAPMAN & HALL/CRC Monographs and Surveys in Pure and Applied Mathematics 
%             R. Luus. Iterative dynamic programming. Chapman and Hall/CRC, 2002.
%
% ************************************************************************

%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'on', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-6),...
 'TolFun', 10^(-6), 'TolCon', 10^(-6), 'MaxFunEval', 15000000, 'MaxIter', 5e10);
% Ode options: see help ode45...

% convergence parameters
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-6, 'AbsTol', 1e-6);

%--------------------------------------------------------------------------
optSPAIDER.conv.f_tol = 0.00047; % tolerance for adaptation algorithm;
optSPAIDER.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
optSPAIDER.conv.t_ref = 256.32; % time for CPU reference;
optSPAIDER.conv.fref = 0.2685; % reference objective functional;
optSPAIDER.conv.tunn = 64; % number of stages
optSPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
optSPAIDER.prob.t0 = 0; % Initial time for integration
optSPAIDER.prob.tf = 5; % Final time for integration
optSPAIDER.prob.x0 = [0 1 0]; %state Initial Condition
optSPAIDER.prob.u_lb = [-1]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [1]; % control vector upper boundary 
optSPAIDER.prob.u0 = [-1]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 1; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF = 7;  % Final Iteration 
%--------------------------------------------------------------------------

%integration options

%--------------------------------------------------------------------------
optSPAIDER.solver.name = 'emso';
optSPAIDER.solver.wfact = 0.4;
optSPAIDER.solver.type = 'external';
optSPAIDER.solver.reg = 0;
optSPAIDER.solver.est = 1;
optSPAIDER.solver.frozen = 0;
optSPAIDER.solver.typ = 0;
optSPAIDER.casestudy.state = 0;
optSPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------

% plot options
%--------------------------------------------------------------------------
optSPAIDER.plot.control = 1;   % 1 for plot or zero for no plot.
optSPAIDER.plot.state = 1; % for plot state
optSPAIDER.plot.plot_after = 0; % plotting storage .mat
optSPAIDER.plot.tf_p = 5; % Final time for integration for plotting
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
optSPAIDER.wavelets.Wav   = 1;  % 1 - Wav; 2 -> ADE...
optSPAIDER.wavelets.adap  = 1; % adaption on(thresholding wavelets) or off (fixed)
optSPAIDER.wavelets.normfrac  = 0.08; % norm
optSPAIDER.wavelets.thr  = 1.0e-5; % fixed
%--------------------------------------------------------------------------

[SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)
