% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrosousa@id.uff.br
% UFF - Universidade Federal Fluminense
% Petrobras

% ************************************************************************
%
%                              #Case study UOTR#
% Problem:     Optimal control of the UOTR model  
% References:  Trica (2022) 
%
%
% ************************************************************************

%load('uotr_ss.mat', 'x0', 'u0', 'par'); % Rodando o modelo no estado estacionário para as condições iniciais

% modelo_u0 = u0;
% modelo_par = par;
% clear u0 par
global PP
PP=blackoil;

global modelo_u0 modelo_par modelo_x0
load('uotr_ctrl_ss_v1.mat','x0','u0','par');
modelo_x0=x0;modelo_x0(end)=eps;
modelo_u0=u0;
modelo_par=par;
clear x0 u0 par

%----------------------------------------------------------------------
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'sqp', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-4),...
 'TolFun', 10^(-4), 'TolCon', 10^(-2), 'MaxFunEval', 15000000, 'MaxIter', 3e10);
% Ode options: see help ode45...
global modelo_dFdx modelo_dFdu
Nx=size(modelo_x0,1);
fhand=@(t,x)uotr_ctrl(t,x,modelo_u0,1,[],modelo_par); %fhand=@(t,x)uotr_ctrl(t,x,[0;0],[u0(1);u0(1)]*1e-3,u0,1,[],par);
fval=fhand(0,modelo_x0);
[modelo_dFdx,~]=numjac(fhand,0,modelo_x0,fval,zeros(Nx,1));

Nu=size(modelo_u0(1),1);
fhand=@(t,u)uotr_ctrl(t,modelo_x0,[u;modelo_u0(2:end)],1,[],modelo_par); %fhand=@(t,u)uotr_ctrl(t,x0,[0;0],u,u0,1,[],par);
fval=fhand(0,modelo_u0(1)); %fval=fhand(0,u0(1)*1e-3);
[modelo_dFdu,~]=numjac(fhand,0,modelo_u0(1),fval,zeros(Nu,1)); %[modelo_dFdu,~]=numjac(fhand,0,u0(1),fval,zeros(Nu,1));

% [n, m] = size(A);
% % Número de blocos diagonais desejados
% num_blocks = 6;  % Queremos 5 blocos de A na diagonal
% 
% % Defina a matriz maior B
% B = zeros(num_blocks * n, num_blocks * m);
% 
% % Preencher a matriz B com blocos de A na diagonal
% for i = 1:num_blocks
%     row_idx = (i-1)*n + (1:n);
%     col_idx = (i-1)*m + (1:m);
%     B(row_idx, col_idx) = A;
% end

% ode parametersmodelo_dFdx
%--------------------------------------------------------------------------
optODE = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10);
% optODE = odeset('RelTol',1e-5,   'AbsTol', 1e-6, 'Jacobian',modelo_dFdx);
% optODE = odeset('RelTol',1e-5,   'AbsTol', 1e-6);
%--------------------------------------------------------------------------
% method of solution
% in_SPAIDER.problem.method = 'single shooting'; 
in_SPAIDER.problem.method = 'single shooting';  %(advisory to be paralell if EMSO)
                          % 'simultaneous'
                          % 'indirect'
in_SPAIDER.problem.type = 'optimal control';    
                        % 'dynamic parameter estimation'
                        % 'multiple dynamic objectives'
                        % 'static optimization'
in_SPAIDER.problem.analitic_derivarive = 'sensitivy analisys';   %(If optimset.GradObj is on and/or optimset.GradConstr is on)
% in_SPAIDER.problem.analitic_derivarive = 'automatic diferentiation';
                                        % 'numeric';
in_SPAIDER.problem.algorithm = 'serial';
%                              'parallel'    % works only for multiple shooting and external routines (as EMSO and Hysys)    
in_SPAIDER.problem.tfree = 0;  % if is free time or fixed time
   in_SPAIDER.problem.interpolation =  'stage';  % interpolation of control profile
%  in_SPAIDER.problem.interpolation =  'spline' ;
%  in_SPAIDER.problem.interpolation =  'linear';
% in_SPAIDER.problem.interpolation =  'nearest';
% in_SPAIDER.problem.interpolation =  'stage';
%  in_SPAIDER.problem.interpolation =  'reg_const';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% convergence options
in_SPAIDER.conv.f_tol = 1e-4; % tolerance for adaptation algorithm;
in_SPAIDER.conv.f_tol2 = 0.8e10;  % second order tolerance for adaptation algorithm;
in_SPAIDER.conv.t_ref = 499; % time for CPU reference;
in_SPAIDER.conv.fref = 0.048055; % reference objective functional;
in_SPAIDER.conv.utarg =0; % u as target functional 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%integration and solver options
in_SPAIDER.solver.type = 'external'; % or internal (for matlab routines)
in_SPAIDER.solver.name = 'ode15s'; % hysys, aspen, ode45, ode23s (and all matlab solvers), dasslc.
in_SPAIDER.solver.est = 1; % option for integration for stages (elements)
in_SPAIDER.solver.reg = 0; % option for regularization function (works only when est = 0!)
in_SPAIDER.solver.speed_factor = 0; % option for trying accelerate solution (may prejudice the quality of solution)
in_SPAIDER.casestudy.status = 'off'; % to computate a case study of various thresholds
in_SPAIDER.casestudy.iter = '16';
in_SPAIDER.solver.sampling = 5;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% spaider plot options
in_SPAIDER.plot.control = 'on';   % 1 for plot or zero for no plot.
in_SPAIDER.plot.state = 0; % for plot state
in_SPAIDER.plot.plot_after = 'off'; % plotting storage .mat
in_SPAIDER.plot.tf_p = 3600; % Final time for integration for plotting
in_SPAIDER.plot.record = 'off'; % Record the generated figures  
in_SPAIDER.plot.video = 'off'; % Record a video of control profile refinement  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% adaptive wavelet options
in_SPAIDER.wavelets.wname = 'db1';  % wavelets type (dbX families: db1, db2...db4, 'mayer')
in_SPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure ('rigrsure', 'heursure', 'sqtwolog', 'gcv', 'bayes', 'cvpshrink')
in_SPAIDER.wavelets.sorh  = 's';  % hard or soft threshold option
in_SPAIDER.wavelets.scal  = 'sln'; % scale procedure  sln, mln or one

in_SPAIDER.wavelets.procedure   = 'on';  % on - adaptation wavelets; 2 -> equidistant refinement...
in_SPAIDER.wavelets.thresholding  = 'automatic';  % norm-based; automatic; fixed.
in_SPAIDER.wavelets.norm_threshold  = 0.08; % norm threshold
in_SPAIDER.wavelets.fix_threshold  = 1.0e-5; % fixed threshold
in_SPAIDER.wavelets.max_refinement = 128; % constraints for number of stages (mudar para max_discret
in_SPAIDER.wavelets.refinement_position = 0.42; % [0.1 0.9] location of refinement in each stage 
in_SPAIDER.wavelets.delt = 1.0e-5*0.1;
%--------------------------------------------------------------------------

try
[out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE,SPAIDER)
catch ME
    
        idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match'); 
        
         if strcmp(idSegLast,'') 
                 [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
         else
                 [out_SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE)
         end
   
end
figure(2)
bar(out_SPAIDER.time_ac(1:out_SPAIDER.iter+2))