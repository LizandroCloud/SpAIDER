% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Problem E#
%                    SemiBacth Ishothermal CSTR Reactor
% References:  Schlegel (2004), PhD Thesis; 
% Srinavan et al. (2003);
%
% ************************************************************************

% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'on', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-8),...
 'TolFun', 10^(-8), 'TolCon', 10^(-8), 'MaxFunEval', 15000000, 'MaxIter', 5e0);

% convergence parameters
% ************************************************************************
optODE = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8);
%--------------------------------------------------------------------------
optSPAIDER.conv.f_tol = 0.0001; % tolerance for adaptation algorithm;
optSPAIDER.conv.f_tol2 = 0.8;  % second order tolerance for adaptation algorithm;
optSPAIDER.conv.t_ref = 2.544e+03; % time for CPU reference;
optSPAIDER.conv.fref = 32.686; % reference objective functional;
optSPAIDER.conv.tunn = 48; % number of stages
optSPAIDER.conv.utarg =0; % u as target functional 
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
optSPAIDER.prob.isF = 4  % Final Iteration 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
optSPAIDER.solver.name = 'ode45'; 
optSPAIDER.solver.wfact = 0.5;
optSPAIDER.solver.type = 'internal';
optSPAIDER.solver.reg = 0;
optSPAIDER.solver.est = 1;
optSPAIDER.solver.frozen = 0;
optSPAIDER.solver.typ = 0;
optSPAIDER.casestudy.state = 0  ;
optSPAIDER.casestudy.iter = '1';
%--------------------------------------------------------------------------


% plot options
%--------------------------------------------------------------------------
optSPAIDER.plot.control = 0;   % 1 for plot or zero for no plot.
optSPAIDER.plot.state = 0; % for plot state
optSPAIDER.plot.plot_after = 0; % plotting storage .mat
optSPAIDER.plot.tf_p = 15; % Final time for integration for plotting
%--------------------------------------------------------------------------


% wavelet options
%--------------------------------------------------------------------------
%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa


% wname = 'db1';  % wavelets type
% tptr = 'heursure'; % thresholding procedure
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
optSPAIDER.wavelets.adap  = 2; % adaption on(thresholding wavelets) or off (fixed)
optSPAIDER.wavelets.normfrac  = 0.45; % norm
optSPAIDER.wavelets.thr  = 4e-1; % fixed
optSPAIDER.wavelets.delt = 0.05; % variation of threshold between iterations
%--------------------------------------------------------------------------

[SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)
