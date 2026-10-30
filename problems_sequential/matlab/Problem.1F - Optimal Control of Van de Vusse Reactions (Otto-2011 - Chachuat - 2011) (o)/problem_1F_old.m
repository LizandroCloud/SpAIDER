% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Case study 1F#
% Problem:     Cinética de Van der Vusse 
% References:  Otto (2010) (text in Portuguese) Msc. Dissertation
% Description: O objetivo de produção neste estudo é maximizar a concentração do componente
%B no tempo final, que é o produto desejado, reduzindo a produção dos componente nãodesejados
%C e D.
% ************************************************************************

% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'on', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-6),...
 'TolFun', 10^(-6), 'TolCon', 10^(-6), 'MaxFunEval', 5e1, 'MaxIter',  3e10);
% matlabpool open 2
% options = optimset('UseParallel','always');
% convergence parameters
% ************************************************************************
optODE = odeset( 'RelTol', 1e-5, 'AbsTol', 1e-5);
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
optSPAIDER.prob.tf = 6; % Final time for integration
optSPAIDER.prob.x0 = [5.5 0 0 0 378]; %state Initial Condition
optSPAIDER.prob.u_lb = [25 300]; % control vector lower boundary 
optSPAIDER.prob.u_ub = [50 700]; % control vector upper boundary 
optSPAIDER.prob.u0 = [25 355]; % scalar initial estimative of control profile for each variable
optSPAIDER.prob.nU = 2; % Number of control variables...
optSPAIDER.prob.is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
optSPAIDER.prob.isF = 5;  % Final Iteration 
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
optSPAIDER.plot.control = 1;   % 1 for plot or zero for no plot.
optSPAIDER.plot.state = 1; % for plot state
optSPAIDER.plot.plot_after = 0; % plotting storage .mat
optSPAIDER.plot.tf_p = 4; % Final time for integration for plotting
%--------------------------------------------------------------------------


% wavelet options
%--------------------------------------------------------------------------
%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa

%--------------------------------------------------------------------------
optSPAIDER.wavelets.wname = 'db1';  % wavelets type
optSPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure
optSPAIDER.wavelets.sorh  = 'h';  % hard or soft threshold option
optSPAIDER.wavelets.scal  = 'sln'; % scale procedure
optSPAIDER.wavelets.Wav   = 1;  % 1 - Wav; 2 -> ADE...
optSPAIDER.wavelets.adap  = 1; % adaption on(thresholding wavelets) or off (fixed)
optSPAIDER.wavelets.normfrac  = 0.45; % norm
optSPAIDER.wavelets.thr  = 4e-1; % fixed
optSPAIDER.wavelets.delt = 0.05; % variation of threshold between iterations
%--------------------------------------------------------------------------

[SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)
