% ************************************************************************
% DSc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
%                              #Problem B#
%                    The mixed catalyst problem
% References:   (Bell and Sargent. 2000) 
% 
%
% ************************************************************************
% Optimization options: see help fmincon...
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'on', 'GradConstr', 'on',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-12),...
 'TolFun', 10^(-12), 'TolCon', 10^(-12), 'MaxFunEval', 15000000, 'MaxIter', 5e0);
% Ode options: see help ode45...
optODE = odeset( 'RelTol', 1e-11, 'AbsTol', 1e-11);

% convergence parameters
% ************************************************************************
f_tol = 1e-6; %tolerance for adaptation algorithm;
f_tol2 = 0.01;
t_ref = 4051; % time for CPU reference;
fref=0.0480554;
fref= 0.0207312;
tunn = 64; % tunning parameter


%--------------------------------------------------------------------------
optSPAIDER.conv.f_tol = 0.0001;
optSPAIDER.conv.f_tol2 = 0.8;
optSPAIDER.conv.t_ref = 1000;
optSPAIDER.conv.fref = 0.431721;
optSPAIDER.conv.tunn = 64;
optSPAIDER.conv.utarg =1;
%--------------------------------------------------------------------------



% ************************************************************************
%problem parameters
t0 = 0; % Initial time for integration
tf = 1; % Final time for integration
x0 =  [1 0];  %state Initial Condition
u_lb = [0]; % control vector  boundary 
u_ub = [1]; % control vector upper boundary 
u0 = [0.5]; % scalar initial estimative olowerf control profile for each variable
nU = 1; % Number of control variables...
is0 = 3;  % Initial Iteration (Number of stages is ns^is0 )...
isF = 5;  % Final Iteration 


%--------------------------------------------------------------------------
optSPAIDER.prob.t0 = 0.0001;
optSPAIDER.prob.tf = 0.0001;
optSPAIDER.prob.x0 = 0.0001;
optSPAIDER.prob.u_lb = 0.0001;
optSPAIDER.prob.u_ub = 0.0001;
optSPAIDER.prob.u0 = 0.0001;
optSPAIDER.prob.nU = 0.0001;
optSPAIDER.prob.is0 = 0.0001;
optSPAIDER.prob.isF = 0.0001;

if optSPAIDER.conv.utarg == 0;
   optSPAIDER.prob.uref=wopt;
   optSPAIDER.prob.tref=tx(:,:,end);
   optSPAIDER.prob.nref=nsx(1,1,end);
else
   optSPAIDER.prob.uref=[1];
    optSPAIDER.prob.tref=[1];
    optSPAIDER.prob.nref=[1];
end
%--------------------------------------------------------------------------


% ************************************************************************
%integration options
frozen=0;
solver='ode45';  % fixed (rungekutta with step (tf-t0)/dp), ode (for ode45, ode23s, ode15 etc. ),  dasslc, hysys or emso.
wfact = 0.5;  % position of new adapt point.
typ=0;  % free or fixed time
reg=0; % regularization
est=1; % stages or full integration

%--------------------------------------------------------------------------
optSPAIDER.solver.name = 'ode45';
optSPAIDER.solver.wfact = 0.3;
optSPAIDER.solver.type = 'internal';
optSPAIDER.solver.reg = 0;
optSPAIDER.solver.est = 1;
optSPAIDER.casestudy.state = 'on';
optSPAIDER.casestudy.iter = '16';
%--------------------------------------------------------------------------


% plot options
% ************************************************************************
tf_p = 1; % Final time for integration for plotting
plot=0;  % 1 for plot or zero for no plot.
plot_state=0; % for plot state
plot_after=0;


%--------------------------------------------------------------------------
optSPAIDER.plot.control = 1;
optSPAIDER.plot.state = 0;
optSPAIDER.plot.plot_after = 0;
optSPAIDER.plot.tf_p = 250;
%--------------------------------------------------------------------------

% ************************************************************************

%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa


% wavelet options
% ************************************************************************
wname = 'db1';  % wavelets type
tptr = 'sqtwolog'; % thresholding procedure
sorh = 'h';  % hard or soft threshold option
scal = 'sln'; % scale procedure
Wav = 1;  % 1 - Wav; 2 -> ADE...
adap=0; % adaption on(thresholding wavelets) or off (fixed)
% ************************************************************************

%--------------------------------------------------------------------------
optSPAIDER.wavelets.wname = 'db1';  % wavelets type
optSPAIDER.wavelets.tptr  = 'heursure'; % thresholding procedure
optSPAIDER.wavelets.sorh  = 's';  % hard or soft threshold option
optSPAIDER.wavelets.scal  = 'sln'; % scale procedure
optSPAIDER.wavelets.Wav   = 2;  % 1 - Wav; 2 -> ADE...
optSPAIDER.wavelets.adap  = 1; % adaption on(thresholding wavelets) or off (fixed)
optSPAIDER.wavelets.normfrac  = 0.1; % fixed
optSPAIDER.wavelets.thr  = 0.6; % fixed
%--------------------------------------------------------------------------

if plot_after==0 
    
for kg = 1:5
disp(sprintf('CASE STUDY: %d', kg));  % writing on the terminal    
if adap==0 % fixed
    tptr1= strcat('fix_heuristic_2_',num2str(optSPAIDER.wavelets.thr));
elseif adap==1
    tptr1= tptr;
elseif adap==2  % norm-based
     tptr1= strcat('norm_based_1_',num2str(optSPAIDER.wavelets.normfrac));
end

if Wav==2
     tptr2=strcat('fixed_',int2str(is0),'_',int2str(isF));
     tptr1= tptr2;
end

file=strcat(tptr1 , '.txt1');
folder = strcat(pwd , '\','results\step_integration_',int2str(est),'\regularization_',int2str(reg),'\',tptr1,'\',file);
Ax=exist(folder,'file');
if Ax~=0
delete(folder);
end
diary(file);
diary on;
% 
%--------------------------------------------------------------------------
[wopt, erropt, grad,hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv,detais0,detais1] = spaider2_(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,optNLP,optODE,tunn,fref,typ,wfact,reg,est,optSPAIDER);
% 
 diary off;

file1=strcat(pwd , '\', file);
folder = strcat(pwd , '\','results\step_integration_',int2str(est),'\regularization_',int2str(reg),'\',tptr1,'\');
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(folder);
[SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
if SUCCESS == 0
    mkdir(folder);
    [SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
end

file2=strcat(folder,tptr1,'.mat');
save (file2) ;

convhist(kg,1:length(conv)) = conv;
thf(kg) = optSPAIDER.wavelets.thr;
timenorm(kg) = time_ac(length(tx(end,end,:)))/t_ref;
optSPAIDER.wavelets.normfrac = thf(kg)  + 1e-1;
optSPAIDER.wavelets.thr = thf(kg)  + 1e-1;
end
else 
plot=1;
for i=2:length((nsx(1,1,:)))
    for k=1:nU
        ncz=(nsx(k,:,i));
        nz = ncz(ncz~=0 ); 
        if isempty(nz)
            nsz=1;
        else
            nsz=nz;
        end
        n_est(k,i-1) = nsx(k,1,i);
         wx(:,:,i)=wx(:,:,i);
    end
end
 [ns] = plotp(t0,tf,is0,i+1,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,wx,tx,n_est,Fobj,time,time_ac,wfact,tunn);

end


delete *.asv;