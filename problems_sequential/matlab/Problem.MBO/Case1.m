% ************************************************************************
% DSc. Lizandro de Sousa Santos 
% Email: lizzandrosousa@id.uff.br
% Universidade Federal Fluminense

% ************************************************************************
%                              #Case study 1.A#
% Problem:SemiBacth Ishothermal CSTR Reactor  
% References:  Srinivasan et al. (2003a) 
%              Schlegel et al. (2005)
%
% ************************************************************************
startup % loadind dafault configurations

%--------------------------------------------------------------------------
% Optimization options: see help fmincon...
%--------------------------------------------------------------------------
optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'on', 'GradConstr', 'on',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 15000000, 'MaxIter', 10e10);

%--------------------------------------------------------------------------
% Matlab ode parameters
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-9, 'AbsTol', 1e-9);

%--------------------------------------------------------------------------
% % Dynamic optimization configurations
%--------------------------------------------------------------------------
in_Crodys.problem.method = 'single shooting'; 
                          % 'multiple shooting'  (advisory to be paralell if EMSO)  
                          % 'simultaneous'
                          % 'indirect'
in_Crodys.problem.type = 'optimal control';    
                        % 'dynamic parameter estimation'
                       
                        % 'static optimization'
in_Crodys.problem.analitic_derivarive = 'sensitivy analisys';   %(If optimset.GradObj is on and/or optimset.GradConstr is on)
% in_Crodys.problem.analitic_derivarive = 'automatic diferentiation';
                                        % 'numeric';
in_Crodys.problem.algorithm = 'serial';
%                              'parallel'    % suggered for only for multiple shooting and external routines (as EMSO and Hysys)    
in_Crodys.problem.tfree = 0;  % if is free time or fixed time
   in_Crodys.problem.interpolation =  'stage';  % interpolation of control profile
%   in_Crodys.problem.interpolation =  'spline' ;
%    in_Crodys.problem.interpolation =  'linear';
% in_Crodys.problem.interpolation =  'nearest';
% in_Crodys.problem.interpolation =  'stage';
%   in_Crodys.problem.interpolation =  'reg_const';
                                  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_Crodys.conv.f_tol = 1e-5; % tolerance for adaptation algorithm;
in_Crodys.conv.f_tol2 = 10;  % second order tolerance for adaptation algorithm;
in_Crodys.conv.t_ref = 3072.36967484694; % time for CPU reference;
in_Crodys.conv.fref = 0.431721; % reference objective functional;
in_Crodys.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_Crodys.solver.type = 'internal'; % or internal (for matlab routines)
in_Crodys.solver.name = 'ode45'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_Crodys.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_Crodys.solver.sampling = 0; % number of sampling points at each stage
in_Crodys.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_Crodys.casestudy.status = 'off'; % to computate a case study of various thresholds
in_Crodys.casestudy.iter = '16';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Crodys plot options
in_Crodys.plot.control = 'on';   % 1 for plot or zero for no plot.
in_Crodys.plot.state = 0; % for plot state
in_Crodys.plot.plot_after = 'off'; % plotting storage .mat
in_Crodys.plot.tf_p = 250; % Final time for integration for plotting
in_Crodys.plot.record = 'off'; % Record the generated figures  
in_Crodys.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_Crodys.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_Crodys.wavelets.tptr  = 'cvpshrink'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_Crodys.wavelets.sorh  = 'h';  % hard or soft threshold option
in_Crodys.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_Crodys.wavelets.procedure   = 'on';  % on - adaptation wavelets; 2 -> equidistant refinement...
in_Crodys.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_Crodys.wavelets.norm_threshold  = 0.8; % norm threshold
in_Crodys.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_Crodys.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_Crodys.wavelets.refinement_position = 0.44; % [0.1 0.9] location of refinement in each stage 
in_Crodys.wavelets.delt = 1.0e-55*0.1;
in_Crodys.wavelets.rel = 'on';
%--------------------------------------------------------------------------

try
[out_Crodys] = spaider_main(in_Crodys,optNLP,optODE,Crodys)
catch ME
    
        idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match'); 
        
         if strcmp(idSegLast,'') 
                 [out_Crodys] = spaider_main(in_Crodys,optNLP,optODE)
         else
                  [out_Crodys] = spaider_main(in_Crodys,optNLP,optODE)
         end
   
end
