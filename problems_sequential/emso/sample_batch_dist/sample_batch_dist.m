% ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Problem EMSO1#
%                    sample_batch_dist.mso
%
% ************************************************************************

% ************************************************************************
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-9),...
 'TolFun', 10^(-9), 'TolCon', 10^(-9), 'MaxFunEval', 15000000, 'MaxIter', 5e10);
% Ode options: see help ode45...

% convergence parameters
% ************************************************************************
optODE = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8);
f_tol = 0.001; %tolerance for adaptation algorithm;
f_tol2 = 0.8;
t_ref = 1000; % time for CPU reference;
fref=0.21;
tunn = 16;

% ************************************************************************
%problem parameters
t0 = 0; % Initial time for integration
tf = 350; % Final time for integration
%tf_p = 15; % Final time for integration for plotting
x0 = [ 298 2.5 0.2 ]; %state Initial Condition
u_lb = [  3.5e6 100]; % control vector lower boundary 
u_ub = [ 4.8e6 1000]; % control vector upper boundary 
u0 = [  3.7743e6 500]; % scalar initial estimative of control profile for each variable
nU = 1; % Number of control variables...
is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
isF = 3;  % Final Iteration 

% ************************************************************************
%integration options
frozen=0;
optSPAIDER.solver.name = 'emso';
solver='emso';  % ,fix_runge, ode (for ode45, ode23s, ode15 etc. ),  dasslc, hysys or emso.
wfact = 0.5; % position of new adapt point.
typ=0;  % free or fixed time
reg=0; % regularization
est=1; % stages or full integration

% plot options
% ************************************************************************
tf_p = 1; % Final time for integration for plotting
plot=1;  % 1 for plot or zero for no plot.
plot_state=1; % for plot state

% ************************************************************************


% ************************************************************************

%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa

% wavelet options
% ************************************************************************
wname = 'db1';  % wavelets type
tptr = 'sqtwolog'; % thresholding procedure
sorh = 's';  % hard or soft threshold option
scal = 'mln'; % scale procedure
Wav = 2;  % 1 - Wav; 2 -> ADE...
adap=1; % adaption on(thresholding wavelets) or off (fixed)

file=strcat(tptr , '.txt1');
folder = strcat(pwd , '\','results\step_integration_',int2str(est),'\regularization_',int2str(reg),'\',tptr,'\',file);
Ax=exist(folder,'file');
if Ax~=0
delete(folder);
end
diary(file);
diary on;
optCONV=1;optWAV=1;optPLOT=1;
 [wopt, erropt, grad,hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv] = spaider3_(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,...
 tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,...
 f_tol2,t_ref,optNLP,optODE,tunn,fref,typ,wfact,reg,est,optCONV,optWAV,optSPAIDER,optPLOT);

diary off;

diary off;

file1=strcat(pwd , '\', file);
folder = strcat(pwd , '\','results\step_integration_',int2str(est),'\regularization_',int2str(reg),'\',tptr,'\');
% [SUCCESS,MESSAGE,MESSAGEID] = mkdir(folder);
[SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
if SUCCESS == 0
    mkdir(folder);
    [SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
end
% for i=2:length((nsx(1,1,:)))
%     for k=1:nU
%         ncz=(nsx(k,:,i));
%         nz = ncz(ncz~=0 ); 
%         if isempty(nz)
%             nsz=1;
%         else
%             nsz=nz;
%         end
%         n_est(k,i-1) = nsx(k,1,i);
%          wx(:,:,i)=wx(:,:,i);
%     end
% end
%  [ns] = plotp(t0,tf,is0,i+1,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,wx,tx,n_est,Fobj,time,time_ac,wfact);

