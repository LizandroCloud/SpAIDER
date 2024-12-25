function[wopt, tssx, erropt, out_jac,out_hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv,detais00,...
    detais11,uth,ld0,ld1,md00,md11,vd00,vd11,difU,difd,mseU,is,condh,tst,fst,out_SPAIDER] =  spaider(in_SPAIDER,optNLP,optODE)

% [out_SPAIDER.control
%  out_SPAIDER.states
%  out_SPAIDER.optimization
%  out_SPAIDER.convergence
%  out_SPAIDER.wavelets
%] =  spaider(in_SPAIDER,optNLP,optODE)


 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrosousa@id.uff.br
 % home page: www.ceopt.com
 % Universidade Federal Fluminense - Departamento de Engenharia Quimica -
 % Niteroi - Brazil
 % ************************************************************************
 %                              Description
 % This code is applyed to solve dynamic optimization problems with
 % wavelets basis as adaptive strategy for mesh discretization. In 
 % this strategy wavelets are used to identify prospective locations
 % at control vector parameterization where de optimization algorithm
 % will refine the discretization. To do so, details coefficents of
 % wavelets are thresholded according to a specific metric, and the
 % thesholded coefficients indicates these locations. For more information,
 % please read the article:
 % http://www.sciencedirect.com/science/article/pii/B9780444595201501032
 % ************************************************************************
 
 % changes:
 % 07-27-2012: organization and generalisation for nU control variables...
 % 08-14-2012: incorporation of tx and w...
 % 09-20-2012: improving meshref function: consideration of bounds effects
 % in refinement...
 %             improving the wavegraf function: incorporation of
 %             approximations and details pictures...
 % 11-25-2012: incorporporation of dyn.m and scipt.m for automatic derivation. 
 % 01-18-2013: incorporation of write_obj.m and write_cr.m
 % 01-19-2013: 
 % 04-10-2013: elimination of randomic effect in refinement.
 % 06-13-2013: storage of details at each iteration
 % 07-20-2013: algebraic vriables printing
 % 09-30-2013: implementing multiples shlooting
 % 08-29-2014: initialization and parameter estmation
 %*************************************************************************
 [wopt, erropt, out_jac,out_hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv,detais00,...
    detais11,uth,ld0,ld1,md00,md11,vd00,vd11,difU,difd,mseU,is,condh,tst, fst] = initialization;


 % first of all, calling the input_model
 [ dx y par e_var a_var c_var J cr_eq cr_deq in_SPAIDER.input] = input_model; % routine where are written the model and optimization
 in_SPAIDER.equations = [dx y par e_var a_var c_var J cr_eq cr_deq];  % computing the equations (in strings format)
 
 [t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,optSPAIDER.wavelets.procedure,adap,frozen,...
    solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,tunn,fref,typ,wfact,reg,est,folder] = spaider_head(in_SPAIDER);
 scal0=scal;
 tssx=[];
 w01=[]; % pre-allocating
 w02=[]; detais00=[]; detais11=[]; uth=[]; ld0=[]; ld1=[]; md00=[]; md11=[]; difU=[]; difd=[]; mseU=[]; vd00=[]; vd11=[]; ts=[];
 tst=1; r0 =1;  fst=1; u_opt=1; li=0;  ci=0;  w0ns=[];  j =0;  sc=0;
 ns=ones(nU,1);
 wL=ones(nU,1);
 wU=ones(nU,1);
 w0=ones(nU,1);
 wLL=ones(1,nU*1);
 wUU=ones(1,nU*1);
 wopti=ones(nU,1);
 wxi=ones(nU,1);
 w=ones(nU,1);
 tx=ones(1,1,1);
 ts=ones(1,1,1);
 tss=ones(1,nU*1);
 f=[];
 thrs = []; % matrix of calculated thresholds 
 kit = 1;
 iter=[1:40];
 Fobj=zeros(1,40);
 time=zeros(1,40);
 time_ac=zeros(1,40);
 MSE=zeros(1,40);
 SNR=zeros(1,40);
 toc_ant=0; 
 conv(1) =40;
 cvc1=[0];
 cdc1=[0];
 convv = 1;
 d_cr_eq=[]
 d_cr_deq=[]
 param=[]
 

 % checking initial configuration...
 
 if in_SPAIDER.plot.tf_p ~= in_SPAIDER.input.tf
                                %Advice
        %******************************************************************
        adv{1}= ' The time plotting is different from the final time integration.';
        adv{2}= ' You should setting the values to be equal.';
        status = 'warning';
        [mess] = spaider_message(adv,status);
        if mess == 1
            return;
        else
            
        end
        %******************************************************************
        end
        
 
 if exist('alg_eq.m') % checking if exists algebraic equations function file .m
 else
  [s]=write_alg('', '', '', 1, 1, 1, 0, ''); % algebraic equations
 end

 if strcmp(in_SPAIDER.solver.name,'emso') && strcmp(in_SPAIDER.solver.type,'internal')  % for emso must be external
        
        %Advice
        %******************************************************************
        adv{1}= ' If the choosen solver is EMSO, the solver.type must be external.';
        adv{2}= ' Changing automatically to external.';
        status = 'advice';
        [mess] = spaider_message(adv,status);
        if mess == 1
            return;
        else
            
        end
        %******************************************************************
        in_SPAIDER.solver.type = 'external'; % changing to external
 else
%      [];
 end

 if strcmp(in_SPAIDER.solver.type,'external')
%        [s]=write_alg('', '', '', 1, 1, 1, 0, ''); % algebraic equations
 else
 % dynamic functions
 % ************************************************************************
 [param eq eq_sens n_p n_c n_a y J dJ cr_eq d_cr_eq cr_deq d_cr_deq cr_eq_type cr_deq_type n_s] = dyn; % calculating Jacobian
 
        
        % checking initial conditions
        
        if length(in_SPAIDER.input.x0) ~= n_s
                                %Error
        %******************************************************************
        adv{1}= ' The number of initial conditions is inconsistent';
        adv{2}= ' Interrupting algorithm...';
        status = 'error';
        [mess] = spaider_message(adv,status);
        if mess == 1
            return;
        else
            
        end
        %******************************************************************
        end
 
 
 
 % writing program routines that depens of Jacobian calculation
 if strcmp(in_SPAIDER.problem.method,'single shooting')
    [s]=write_alg(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % algebraic models
    [s]=write_models(in_SPAIDER,eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % dynamic models
    [s]=write_obj(J, dJ,param); % objective
    [s]=write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints
 elseif strcmp(in_SPAIDER.problem.method,'multiple shooting')
    [s]=write_alg(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % algebraic models
    [s]=write_models(in_SPAIDER,eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % dynamic models
    [s]=write_obj(J, dJ,param); % objective
    [s]=write_cr_ms(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints
    % [s]=write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints

 end
 
 %*************************************************************************
  end
 
  
 %            Detecting wavelets algorithm
 %*************************************************************************
 if strcmp(in_SPAIDER.wavelets.procedure,'on')   % if uses wavelets
    if strcmp(in_SPAIDER.wavelets.thresholding,'automatic')  % if wavelet automatic threshold 
        file = in_SPAIDER.wavelets.tptr; 
    elseif strcmp(in_SPAIDER.wavelets.thresholding,'fixed') % if wavelet fixed threshold 
        file = 'fixed threshold';
    elseif strcmp(in_SPAIDER.wavelets.thresholding,'norm-based') % if wavelet norm threshold 
        file = 'norm threshold';
    end
 else % if no wavelets 
     file = 'no-wavelets'; 
 end
 %*************************************************************************

% for mulptible shooting
if strcmp(in_SPAIDER.problem.method,'multiple shooting')
            
    if isempty(in_SPAIDER.input.x0_lb) && length(in_SPAIDER.input.x0_lb)< length(in_SPAIDER.input.x0)
                                                                   
        %******************************************************************
        adv{1}= ' For multiple shooting you must specify the constraints of initial values';
        adv{2}= ' Please, go to routine: model_input.m and specify data.x0_lb and data.x0_ub';
        status = 'error';
        [mess] = spaider_message(adv,status);
        if mess == 1
            return;
        else
            
        end
    end
    
    
    if typ==1  % if time free
    nUt = nU;  % original number of control variables
    nU = nUt + length(in_SPAIDER.input.x0) + 1;  % number of controls + initials + time
    else
    nUt = nU;  % number of controls 
    nU = nUt + length(in_SPAIDER.input.x0);  % number of controls + initials 
    end
else
    if typ==1  % if time free
    nUt = nU;  % original number of control variables
    nU = nUt + 1;  % number of controls +  time
    else
    nUt = nU;
    end
end
    
% if typ == 1   % flag that specifies if is time free or not
%     nUt = nU-1;
% else
%     nUt=nU;
% end 
 
for is = is0:isF  % number of stages is ns^is0  % main loop...
   
    if is==is0   % Initial (first guess) configuration...
        for j = 1:nUt; % for each control, variable
        ns(j) = (2^is); % initial number of stages
        wL(j,1:ns(j)) = u_lb(j)*ones(ns(j),1);  % lower boundary of u
        wU(j,1:ns(j)) = u_ub(j)*ones(ns(j),1);  % upper boundary of u
        w0(j,1:ns(j)) = u0(j)*ones(ns(j),1);  % initial estimative da u...
        ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
        w_frozen = w0;
        check = [2 3 4 5 6 7 8 9];
        ts_frozen = ts;
        ts_a = ts;
        end
%     if strcmp(in_SPAIDER.problem.method,'multiple shooting')
%         for j = 1+nUt:nUt+nU-1; % for each control, variable
%         wL(j,1:1) = in_SPAIDER.input.x0_lb(j-nUt)*ones(1,1);  % lower boundary of u
%         wU(j,1:1) = in_SPAIDER.input.x0_ub(j-nUt)*ones(1,1);  % upper boundary of u
%         w0(j,1:1) = in_SPAIDER.input.x0(j-nUt)*ones(1,1);  % initial estimative da u...
%         end     
%     end
    
        
        
    end
    ki=1; k=0;
%     w00=ones(1,nU*1);
    tss=[];
    ns_a = ns(1:nUt);
    for i=1:nUt
         for j=k+1:ns_a(i)+k
             tss(j)=ts(i,ki);
             ki=ki+1;
         end
         tss(j+1)=ts(i,ki);
         k = j+1;
         ki=1;
    end
   
    ts_ms=unique(tss);
    if strcmp(in_SPAIDER.problem.method,'multiple shooting')
        for j = 1+nUt:nU; % for each control, variable
        wL(j,1:length(ts_ms)-1) = in_SPAIDER.input.x0_lb(j-nUt)*ones(length(ts_ms)-1,1);  % lower boundary of u
        wU(j,1:1:length(ts_ms)-1) = in_SPAIDER.input.x0_ub(j-nUt)*ones(length(ts_ms)-1,1);  % upper boundary of u
%         w0(j,1:1:length(ts_ms)-1) = in_SPAIDER.input.x0(j-nUt*1,1:length(ts_ms)-1)*ones(length(ts_ms)-1,1);  % initial estimative da u...
        if is==is0
             % w0(j,1:1:length(ts_ms)-2) = in_SPAIDER.input.x0(1,1:length(ts_ms)-2);
             w0(j,1:1:length(ts_ms)-1) = in_SPAIDER.input.x0(1,j-nUt)*ones(length(ts_ms)-1,1);
              
        else
            w0(j,1:1:length(ts_ms)-2) = in_SPAIDER.input.x0(j-nUt*1,1:length(ts_ms)-2);
        end
        end     
    end
%  iu=is+2; % graphical counter....
    ns_a = ns(1:nUt);
    if typ == 1
            wL(nU,1) =  u_lb(nU);
            wU(nU,1) =  u_ub(nU);
            if is==is0
            w0(nU,1) = u0(nU);  % initial estimative da u...\
            else
            w0(nU,1) =  talf;
            end         
    end

    disp(sprintf('Iteration: %d', is));  % writing on the terminal
    disp(sprintf('Number of stages: %d', ns));
    disp(sprintf('Method: %s',file)) 
    % calling optimization algorithm...
     % pre-optimization#####################################################
    k=0;
     tss(1) = ts(1,1);
% 	ki=1;
%     w00=ones(1,nU*1);
%     tss=[];
%     for i=1:nUt
%          for j=k+1:ns_a(i)+k
%              tss(j)=ts(i,ki);
%              ki=ki+1;
%          end
%          tss(j+1)=ts(i,ki);
%          k = j+1;
%          ki=1;
%     end
%     if strcmp(in_SPAIDER.problem.method,'single shooting')
    k=0; ki=1;
    for i=1:nUt
     for j=k+1:ns(i)+k
         w00n(j)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
        w00(j)=w0(i,ki);
        w00a(j) = w0(i,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
   k=0; ki=1;
   for i=1:nUt
     for j=k+1:ns(i)+k
         wLL(j)=wL(i,ki);
         wUU(j)=wU(i,ki)';
          wLLn(j)=0;
          wUUn(j)=1;
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
   
%     end
%  ko=k;
%  ts_a = ts_frozen;
%  k=0;
%  woi = length(w00);
 if strcmp(in_SPAIDER.problem.method,'multiple shooting')
     ko=k;
 ts_a = ts_frozen;
 k=0;
 woi = length(ts_ms)-1;
    for i=1:nU-1
     for j=k+nUt:length(ts_ms)-3+k+nUt
         w00n(j+woi)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
         w00(j+woi)=w0(i+1,ki);
         w00a(j+woi) = w0(i+1,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
    end
%    woi = length(w00);
    k=0; ki=1;
   for i=1:nU-1
     for j=k+nUt:length(ts_ms)-3+k+nUt
%          w00n(j+woi)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
        wLL(j+woi)=wL(i+1,ki);
        wUU(j+woi) = wU(i+1,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
     
     
%         k=0;
%         j=1:nU*length(ts_ms)-2
%         w00(j+ko) = w0(j+1,1);
%         wLL(j+ko) = wL(j+1,1);
%         wUU(j+ko) = wU(j+1,1);   


   
 end
 
 if typ==1
 nss = sum(ns);
 w0a = w0(nU,1);
%  w00(nss+1) = w0(nU,1)-u_lb(nU)/(u_ub(nU)-u_lb(nU));
 w00(nss+1) = w0(nU,1);
 wLL(j+1) = wL(nU,1);
 wUU(j+1) = wU(nU,1);
%  wLL(j+1)=0;
%  wUU(j+1)=1;
 else
 nss = sum(ns);
 end

%  if is==5
% w0(2,:)=[ 
% 2.0672E+01
% 2.0168E+01
% 1.0084E+02
% 1.0084E+02
% 1.0084E+02
% 8.9748E+01
% 7.9664E+01
% 7.0588E+01
% 6.2017E+01
% 5.3950E+01
% 4.7395E+01
% 4.1849E+01
% 3.6303E+01
% 3.3782E+01
% 1.0034E+02
% 1.0034E+02
% 1.0034E+02
% 1.0034E+02
% 1.0034E+02
% 2.0672E+01
% 2.0672E+01
% 2.1176E+01
% 2.1176E+01
% 2.1176E+01
% 2.1176E+01
% 2.0672E+01
% 1.0034E+02
% 1.0034E+02
% 9.9832E+01
% 9.9832E+01
% 1.0034E+02
% 9.9832E+01
% 
% ];
% 
% w0(1,:)=[ 
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% ];
%  end


 % optimization
 %########################################################################
%  [f_inicial] = objfunc(x0,nss,tss,w00a,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,in_SPAIDER)

%  objective = @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
%  nonl_constraint = @(ws)constr(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
%  if is==2
%       in_SPAIDER.solver.sampling=9;
%  elseif is==3
%       in_SPAIDER.solver.sampling=4;
%  elseif is==5
%       in_SPAIDER.solver.sampling=1;
%  end

 switch (in_SPAIDER.solver.name)  % for each solver
     case 'mbo'
     in_SPAIDER.solver.func = @ode15i;
     objective = @(ws)objfunc_mbo(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_mbo(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'unisim'
     in_SPAIDER.solver.func = @unisim;
     objective = @(ws)objfunc_unisim(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_unisim(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'emso'
     in_SPAIDER.solver.func = @emso;
     objective = @(ws)objfunc_emso(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_emso(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode45'   
      in_SPAIDER.solver.func = @ode45;
      objective = @(ws)objfunc_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
      case 'ode15i'   
      in_SPAIDER.solver.func = @ode15i;
      objective = @(ws)objfunc_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode23s'   
     in_SPAIDER.solver.func = @ode23s;    
     objective = @(ws)objfunc_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode15s'
     in_SPAIDER.solver.func = @ode15s;      
     objective = @(ws)objfunc_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode23t'
     in_SPAIDER.solver.func = @ode23t;       
     objective = @(ws)objfunc_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'dasslc' 
     in_SPAIDER.solver.func = @dasslc;
     objective = @(ws)objfunc_dasslc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_dasslc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER); 
 end
 % calling optimization routine...
%  
%   if strcmp(in_SPAIDER.problem.method,'multiple shooting')
%       if is==is0
%             optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
%  'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
%  'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 15000000, 'MaxIter', 5);
%       end
%   end
% if is>is0-1
%      in_SPAIDER.solver.sampling = 100;
     
    if strcmp(in_SPAIDER.problem.method,'multiple shooting')  % initical guess for multiple shooting ...
         in_SPAIDER.problem.method='single shooting';
            [r0 tst fst] = g_state(plot_state,nss,tss,w0(1:nUt,:),ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,1,reg,est,wLL, wUU, in_SPAIDER);
         in_SPAIDER.problem.method='multiple shooting';
                % [s]=write_cr_ms(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints

    delta = (length(tst)/(ns));
        x_ms = [];
        for jk=1:length( fst(1,:) )
            xz = []; yz = [];
            xz = tst(:);
            yz = fst(:,jk);
            yi = interp1(xz,yz,'spline', 'pp');  % interpolation
            for jm=1:(ns)-2
                x_ms(jk,jm) = ppval(yi,delta*jm);
            end
        end
%         x_ms(jk,jm) = fst(delta*jm,jk);   
        for jk=1:length(fst(1,:))
            x_ms(jk,(ns)-1) = fst(end,jk);
        end
        
        w0(2:end,2:end)=x_ms;
     in_SPAIDER.input.x0 = x_ms;    
     
      if strcmp(in_SPAIDER.problem.method,'multiple shooting')
     ko=k;
 ts_a = ts_frozen;
 k=0;
 woi = length(ts_ms)-1;
    for i=1:nU-1
     for j=k+nUt:length(ts_ms)-3+k+nUt
         w00n(j+woi)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
        w00(j+woi)=w0(i+1,ki);
        w00a(j+woi) = w0(i+1,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
    end
%    woi = length(w00);
    k=0; ki=1;
   for i=1:nU-1
     for j=k+nUt:length(ts_ms)-3+k+nUt
%          w00n(j+woi)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
        wLL(j+woi)=wL(i+1,ki);
        wUU(j+woi) = wU(i+1,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
     
     
%         k=0;
%         j=1:nU*length(ts_ms)-2
%         w00(j+ko) = w0(j+1,1);
%         wLL(j+ko) = wL(j+1,1);
%         wUU(j+ko) = wU(j+1,1);   
 end
     
     
    end
  
     
     
     
   % end 
% 
% if is==5
% w0(2,:)=[ 
% 2.0672E+01
% 2.0168E+01
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 8.9748E+01
% 7.9664E+01
% 7.0588E+01
% 6.2017E+01
% 5.3950E+01
% 4.7395E+01
% 4.1849E+01
% 3.6303E+01
% 3.3782E+01
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 2.0672E+01
% 2.0672E+01
% 2.1176E+01
% 2.1176E+01
% 2.1176E+01
% 2.1176E+01
% 2.0672E+01
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 1.00E+02
% 
% ];
% 
% w0(1,:)=[ 
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 5.78
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% ];
% w00=[w0(1,:) w0(2,:)];
% end


% wUU=[w0(2,:) w0(1,:)]*1.2;
% wLL=[w0(2,:) w0(1,:)]*0.8;
 tic
 [wopt, erropt, iout,output,lambda,grad,hessian] = ....
     runobjconstr( w00, wLL, wUU, nonl_constraint , optNLP,...
     x0,nss,tss,ns,length(x0),nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,in_SPAIDER);

 % nvars = ns*nUt;

 % [wopt, erropt, iout,output,lambda,grad,hessian] = particleswarm(runobjconstr( w00, wLL, wUU, nonl_constraint , optNLP,...
 %     x0,nss,tss,ns,length(x0),nUt,w_frozen,check,...
 %         ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,in_SPAIDER),nvars,wLL,wUU);

 toc
 tssx = tss;
 
% 
%  tic
%  [wopt, erropt, iout,output,lambda,grad,hessian] = ....
%      fmincon(  objective, w00, [], [], [], [], wLL, wUU, nonl_constraint , optNLP);  
%  toc
%  options = gaoptimset('MutationFcn',@mutationadaptfeasible, 'PlotFcns',...
%          {@gaplotbestf,@gaplotbestindiv,@gaplotexpectation,@gaplotstopping}, 'Generations', 1000,...
%          'TolFun', 10^(-12), 'TolCon', 10^(-12));
%   tic
%  [wopt, erropt, iout] = ga(  objective, ns, [], [], [], [], wLL, wUU, nonl_constraint , options);  
%  toc
%  
 % time counting
 if is==0
     is = is+ 1;
%      break;
 end
 iter(is)=is;
 time(is)=toc;
 time_ac(is)=time(is)+toc_ant;
 toc_ant=time_ac(is);
 kit=kit+1;
 
 % storaging hessian...
 szh=size(hessian);
 condh(is) = cond(hessian); 
 szg=length(grad);
 out_hessian(is,1:szh(1),1:szh(1)) = hessian;
 out_jac(is,1:szg) = grad;
 
 

 
% % correcting normalization
%   k=0; ki=1;
%   for i=1:nUt
%      for j=k+1:ns(i)+k
%          wopt(j)=wopt(j)*(wUU(j)-wLL(j))+wLL(j);
%          ki=ki+1;
%      end
%      k = j;
%      ki=1;
%   end


% [f_final] = objfunc(x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,in_SPAIDER)
 
 % pos optimization********************************************************
 
 if typ==1
    talf=wopt(end);
 else
    talf=1;
 end

% State graphs
if strcmp(in_SPAIDER.problem.type,'parameter estimation')
else
if is>is0-1
    if plot_state==1  % plotting only in this case
     [r0 tst fst] = g_state(plot_state,nss,tss,wopt,ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,talf,reg,est,wLL, wUU,in_SPAIDER);
%      [r0 tst fst] = g_constraints(plot_state,nss,tss,wopt,ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,talf,reg,est,wLL, wUU,in_SPAIDER);
    end
end        
end
     if is==is0
    % guess for multiple shooting
    if strcmp(in_SPAIDER.problem.method,'multiple shooting')
        if (-1)^(length(tst)) == 1
        delta = ((length(tst))/(ns*2));
        elseif (-1)^(length(tst)) == -1
        delta = ((length(tst)-1)/(ns*2));
        end
        x_ms = [];
        for jk=1:length( fst(1,:) )
            xz = []; yz = [];
            xz = tst(:);
            yz = fst(:,jk);
            yi = interp1(xz,yz,'spline', 'pp');  % interpolation
            for jm=1:(ns*2)-2
                x_ms(jk,jm) = ppval(yi,delta*jm);
            end
        end
        x_ms(jk,jm) = fst(delta*jm,jk);   
        for jk=1:length(fst(1,:))
            x_ms(jk,(ns*2)-1) = fst(end,jk);
        end
     in_SPAIDER.input.x0 = x_ms;    
    end
  
 
    end
 
% pause
% %############ reconfigurating wopt for frozen option 
    if frozen==1
    kq=1;
    for i=1:nUt
            j=1:ns(i);
            wopti(i,j)=wopt(kq);
            kq=kq+1;
    end
    
    ko=0;
    if is == is0
     for i=1:nUt
        for j=1:ns_a(i)
                  wxi(i,j)=wopt(ko+1);
                  ko=ko+1;
        end    
     end
    else
    for i=1:nUt
            for j=1:ns_a(i)
                    for k = 1:length(check(i,:))
                        if j == check(i,k)
                            wxi(i,j)=wopti(i,k);
                            break;
                        else
                            wxi(i,j)=w_frozen(i,j);
                        end
                    end
             end
    end
    end
    ts=ts_a;
    wopt=wxi;
    ns=ns_a;
    else
    kq=1;
    for i=1:nUt
        for j=1:ns(i)
            wxi(i,j)=wopt(kq);
            kq=kq+1;
        end
    end
    wxi_ms = wopt; % for multiple shooting
    wopt=[];
    wopt=wxi;
    ts_a=ts;
    end;
%#########################################################
    
    if strcmp( in_SPAIDER.problem.type, 'parameter estimation')  % If a problem is parameter estimation..stop
              
        %Advice
        %******************************************************************
        adv{1}= ' Parameter estimation done';
        adv{2}= ' The values will be shown bellow';
        status = 'warning';
        [mess] = spaider_message(adv,status);
        parameters = wopt;
        parameters
        hold on;
          [r0 tst fst] = g_state_parameter(plot_state,nss,tss,wopt,ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,talf,reg,est,in_SPAIDER);

        
        if mess == 0
            return;
        else
            
        end
     
    end

    ic = 1:nUt;
    wx(ic,1:length(wopt(ic,:)),is)=wopt(ic,:);  % each line of w is a control variable...
    tx(ic,1:length(ts(ic,:)),is)=ts(ic,:); % time stages
    nsx(ic,1,is)=ns(ic);  % number of stages

 for ic = 1:nUt  % nUt is the number of control variables
% for ic = [1,2,3,4]    
    w0_less = w0;
%     if ic==3 
%        tptr='sqtwolog';
%     end
%     

%      if ic==1
%          tptr = 'sqtwolog'; scal = 'sln'; % scale procedure
%      end
%      if ic==3
%          tptr = 'sqtwolog'; scal = 'sln'; % scale procedure
%      end
if is==8
    aq=1;
end
    if strcmp(optSPAIDER.wavelets.procedure,'on') % Wavelets analysis
        if (is>is0) && (is<=isF)    % if the iteration is superior than first iteration and less then last iteration 
            wopt(ic,1:(ns(ic))) = wopt(ic,1:(ns(ic)))'; %transposing
            [cfd,detais0,maxlev,C,L] = wav_mesh(wopt(ic,1:(ns(ic))),is,ic,wname); % into wavelets domain
            % calculating the wavelet entropy
            
            out_SPAIDER.energy(is) = sum(C.^2);
            out_SPAIDER.shannon(is) = wentropy(C,'shannon');
            out_SPAIDER.e2s(is) = wentropy(C,'shannon') / sum(C.^2);
            
            out_SPAIDER.logenergy(is) = wentropy(C,'log energy');
            out_SPAIDER.e2l(is) = wentropy(C,'log energy') / sum(C.^2);
            
            
            detais00(ic,1:length(C),is-is0)=C;  % details before thresholding
            ld0(ic,is-is0)=length(C);
            [md0,vd0,muci,sci] = normfit(C);
            md00(ic,is-is0)=md0;
            vd00(ic,is-is0)=vd0;
            if var(wopt(ic,:))<1
                scal = 'one';
            else
                scal = scal0;
            end
            
            [Cth, thrs1, detais1, difdd, w01, s, threshold] = wav_thresh(C,L,maxlev,wname,detais0,wopt,...
                tptr,sorh,adap,scal,ic,is,is0,in_SPAIDER.wavelets.norm_threshold,in_SPAIDER,ns); % Wavelets thresholding
            uth(ic,1:length(w01),is-is0)=w01(ic,:);
            detais11(ic,1:length(Cth),is-is0)=Cth; % details after thresholding
            ld1(ic,is-is0)=length(Cth);
            [md1,vd1,muci,sci] = normfit(Cth);
            md11(ic,is-is0)=md1;
            vd11(ic,is-is0)=vd1;
            thrs(ic,1:length(thrs1),is-is0)=thrs1; 
            difU(ic,1:ns(ic),is-is0)= wx(ic,1:ns(ic),is-is0) - uth(ic,1:ns(ic),is-is0);
            difd(ic,1:length(C),is-is0)=difdd(ic,1:length(C),end);
            for i=1:ns(ic)
                mseU(ic,is-is0) = (abs(wx(ic,i,is-is0)) - abs(uth(ic,i,is-is0)))^2/length(wopt(ic,:));
            end
            [Pgrad] = adapt(w01,Cth,L,maxlev,ic,tunn,ts,adap,in_SPAIDER.wavelets.norm_threshold) % prospective points...
%       if ic==4 
%        Pgrad=0;
%       end           
            
             Pgrad=Pgrad(Pgrad~=0)
            PgradX(ic)=norm(Pgrad);
            if ( isempty(PgradX(ic))) 
                 elseif PgradX(ic)==0
                 else
              [Pgrad] = npoints(Pgrad,wopt(ic,:));
            end
            
%             if output.iterations < 10 
%                 Pgrad = downsample(Pgrad,2)
%             end
%           stop barrier
%             PgradX(ic)=norm(Pgrad);
            if ( isempty(PgradX(ic))) 
                stop_flag(ic)=0;
            elseif PgradX(ic)==0
                stop_flag(ic)=0;
            else
                stop_flag(ic)=1;
            end 
            if ic==nUt
                stop_flag_norm=norm(stop_flag);
                if ( isempty( stop_flag_norm )) 
                     avc=1
                      out_SPAIDER.psnr = [];
                    return
                elseif   stop_flag_norm ==0
                    avc=1
                    out_SPAIDER.psnr = [];
                    return
                end
            end
                
            w02(ic,: ) = wopt(ic,: ); % backing to original optimal profile...
            if plot ==1
                ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10+1;
                ix3 = ix*100+2;
                ix4 = ix*1000+3;
                ts_x=[];
                C_x=[];
                % generation of wavelet mesh before thresholding
                [r_ix1a] =  g_ix1a(ix,ic,ns,ns_a,cfd,nUt,maxlev);
                % Optimal dynamic control profile picture
                [r_ix2a u_par] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER);
                % Control profile before thresholding
                u_parx(ic,is,:)=u_par;
                [r_ix3a] = g_ix3a(ix3,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER);
                % Details bars before thresholding
                [r_ix2b] = g_ix2b(ix2,ic,ns_a,talf,C,nUt,L);
                 % Details bars after thresholding
                [r_ix3b] = g_ix3b(ix3,ic,ns_a,talf,C,nUt,L);
                % Thresholded control
                [r_ix4a fit] = g_ix4a(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,Cth,nUt,w01,wopt,is,in_SPAIDER); 
                % into wavelets domain
                [cfdth,detais00,maxlev] = wav_mesh(w01,is,ic,wname);
                % generation of wavelet mesh after thresholding 
                [r_ix1b] =  g_ix1b(ix,ic,ns,ns_a,cfdth,nUt,maxlev);
                % Thresholded details histogram
                [r_ix4b] =  g_ix4b(ix4,ic,Cth,is,nUt,difd,ns_a);
            else
                ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10;
                ix3 = ix*100+1;
                ix4 = ix*1000+2;
                ts_x=[];
                C_x=[];
                [r_ix4c fitt] = g_ix4c(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,Cth,nUt,w01,wopt,is); 
            end
 % preparing for mesh refinement...
        else  % if is=is0
            w02(ic,: )=wopt(ic,: ) ;
            Pgrad = [ 2 3 4 5 6 7 8 ];   % Initial refinement
            ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10+1;
                ix3 = ix*100+2;
                ix4 = ix*1000+3;
                ts_x=[];
                C_x=[];
                % generation of wavelet mesh before thresholding
%                 [r_ix1a] =  g_ix1a(ix,ic,ns,ns_a,cfd,nUt,maxlev);
                % Optimal dynamic control profile picture
                C=[];
                [r_ix2a u_par] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER);
        end
    else % no wavelets
    
%         if is>is0  % for first iteration...
%             for i = 1:length(wopt(ic,: ))
%                 wopt(ic,i) = wopt(ic,i) + rand(1)*1e-4;  % forcing swmall variability...
%                 w02(ic,:)=wopt(ic,: );
%             end
%         else
            w02(ic,1:ns(ic))=wopt(ic,1:ns(ic));
            w0ns(ic,1:ns(ic) )=wopt(ic,1:ns(ic));
%         end
        Pgrad=2:length(wopt(ic,: ));
        convv = 0;
        detais0 = 0;
        detais1 = 0;
        if (is>=is0) && (is<=isF) 
        if plot ==1
                ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10+1;
                ix3 = ix*100+2;
                ix4 = ix*1000+3;
                ts_x=[];
                C_x=[]; C=1;
                % generation of wavelet mesh before thresholding
                %[r_ix1a] =  g_ix1a(ix,ic,ns,ns_a,cfd,nUt,maxlev);
                % Optimal dynamic control profile picture
                [r_ix2a u_par] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER);
                u_parx(ic,is,:)=u_par;
                % Control profile before thresholding
                [r_ix3a] = g_ix3a(ix3,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is,in_SPAIDER);
                % Details bars before thresholding
               % [r_ix2b] = g_ix2b(ix2,ic,ns_a,talf,C,nUt,L);
                 % Details bars after thresholding
                %[r_ix3b] = g_ix3b(ix3,ic,ns_a,talf,C,nUt,L);
                % Thresholded control
                %[r_ix4a fitt] = g_ix4a(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,Cth,nUt,w01,wopt,is); 
                % into wavelets domain
               % [cfdth,detais0,maxlev] = wav_mesh(w01,is,ic,wname);
                % generation of wavelet mesh after thresholding 
               % [r_ix1b] =  g_ix1b(ix,ic,ns,ns_a,cfdth,nUt,maxlev);
                % Thresholded details histogram
               % [r_ix4b] =  g_ix4b(ix4,ic,Cth,is,nUt,difd,ns_a);
                [r_ix4c fit] = g_ix4c(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,1,nUt,w01,wopt,is); 
                ix=0;
                ix4=0;
            else
                ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10;
                ix3 = ix*100+1;
                ix4 = ix*1000+2;
                ts_x=[];
                C_x=[];
                [r_ix4c fitt] = g_ix4c(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,1,nUt,w01,wopt,is); 
        end
        end
        
    end
%     ixx = [ix;ix2;ix3;ix4 ];
%     ixy = ixx(ixx~=0);
% if plot ==1 && is>3 
%     if (ix~=0)
%         ixs = [ix;ix2;ix3;ix4 ];
%     else
%         ixs = [ix2;ix3;ix4 ]; 
%     end
% [pic2xls] = spaider_xls(ixy,is,ic,folder);
% end
% grid refinement  (meshref function)...
%##########################################################################

[ts_ref w_ref ns_ref point] = meshref(Pgrad,ns,ic,ts_a,w02,t0,tf,wfact,u_par,in_SPAIDER);
ts(ic,1:ns_ref+1)=ts_ref;
w0ns(ic,1:ns_ref)=w_ref;
ns(ic)=ns_ref;  
    
    % guess for multiple shooting
        if strcmp(in_SPAIDER.problem.method,'multiple shooting')
        delta = (length(tst)/(ns(ic)));
        x_ms = [];
         for jk=1:length( fst(1,:) )
%             jk=ic;
            xz = []; yz = [];
            xz = tst(:);
            yz = fst(:,jk);
            yi = interp1(xz,yz,'spline', 'pp');  % interpolation
            for jm=1:(ns(ic))-2
                x_ms(jk,jm) = ppval(yi,ts(jm+1));
            end
         end
%         x_ms(jk,jm) = fst(delta*jm,jk);   
         for jk=1:length(fst(1,:))
            x_ms(jk,(ns(ic))-1) = fst(end,jk);
         end
      in_SPAIDER.input.x0 = x_ms;    
    end
 



if (is>is0) && (is<isF) 
    if plot ==1
        % Control profile with new points
        
        [r_ix3c u_opt] = g_ix3c(ix3,ic,ns,ns_a,t0,tf_p,ts,ts_a,talf,C,nUt,wopt,w0ns,point,is);
        
        % Threshoding in bars
       
        [r_ix3d] = g_ix3d(ix3,is,is0,ic,thrs,L,C);
     
    else
        ix3=1;
        [r_ix3e u_opt] = g_ix3e(ix3,ic,ns,ns_a,t0,tf_p,ts,ts_a,talf,C,nUt,wopt,w0ns,point,is);
    end
end
wL(ic,1:ns(ic)) = u_lb(ic)*ones(ns(ic),1);  % lower boundary of u
wU(ic,1:ns(ic)) = u_ub(ic)*ones(ns(ic),1);  % upper boundary of u
wref(ic,1:length(w0(ic,1:length( ns (ic) ) )),is-1)=w02(ic,1:length( ns (ic) ) );
if (is>is0)
    ixx = [ix;ix2;ix3;ix4 ];
    ixy = ixx(ixx~=0);
if plot ==1 && is>3 
    if (ix~=0)
        ixs = [ix;ix2;ix3;ix4 ];
    else
        ixs = [ix2;ix3;ix4 ]; 
    end
%  [pic2xls] = spaider_xls(ixy,is,ic,folder);
end

try
% close all;
movefile ('*.png' , folder);
movefile ('*.emf' , folder);
movefile ('*.eps' , folder);
movefile ('*.fig' , folder);
movefile ('*.xlsx' , folder);
end
end
 end
%  in_SPAIDER.input.x0 = x_ms;
% if plot ==1 && is>3 
%     if (ix~=0)
%         ixs = [ix;ix2;ix3;ix4 ];
%     else
%         ixs = [ix2;ix3;ix4 ]; 
%     end
% [pic2xls] = spaider_xls(ixs,is,folder);
% end
 
Fobj(is) = erropt; % Objective function value
if frozen==1  % frozen some constraints...It is not working well...
    check=[];
    w_frozen = [];
    ts_frozen=[];
    w_froz=[];
    wL_froz=w0ns;
    wU_froz=w0ns;
    t_froz = [];
    ns_a=[];
    ns_a = ns;
    for i=1:nUt
        t_froz2=setxor(ts_a(i,:),ts(i,:));
        t_froz=setdiff(ts(i,:),ts_a(i,:));
        if (length( t_froz ) >  length(ts(i,:)) )
            t_froz = ts(i,:);
        end
            for j=1:length(t_froz)
                check(i,j) = find(t_froz(j) == ts(i,:));
                w_froz(i,j)=w0ns(i,check(i,j));
                wL_froz(i,j)=wL(i,check(i,j));
                wU_froz(i,j)=wU(i,check(i,j));
            end
            ns(i) = length(t_froz);
            t_froz = [];
    end
    w_frozen=w0ns;
    ts_frozen=ts;
    wL = []; wU = []; w0ns = [];
    wL = wL_froz;
    wU = wU_froz;
    w0ns = w_froz;
else
    ns_a = ns; 
    w_frozen=w0ns;
    ts_frozen=ts;
end
 li=0;
 ci=0;
 C=0;
 L=0;
 w01=[];
 w02=[];
 ts_a = [];
 w0 = w0ns;
%  w0(4,1:8) = wopt(4,1:8);
 w0ns=[];
 wLL=[];
 wUU=[];
 ts_r =[];
 
 % Stop Criterium 
 if is>is0
 if strcmp(optSPAIDER.wavelets.procedure,'on')
 
  u_parx(is,1:length(u_par) )=u_par;
 [sc cvc cdc] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,u_opt,cvc1,cdc1,time_ac,nUt,u_parx);
  cvc1=cvc;
  cdc1=cdc;
   conv(is)=cvc1(end);
   convv(is)=cdc1(end);
%    if sc==1
%        return
%    end
 
 else
     u_opt=1;
   [sc cvc cdc] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,u_opt,cvc1,cdc1,time_ac,nUt,u_parx);
   cvc1=cvc;
   cdc1=cdc;
   conv(is)=cvc1(end);
   convv(is)=cdc1(end);
 end
  end
%% Plotting Objective Function
if plot==1
    if is>is0+1
% figure(1) % plot objective function at each iteration
% subplot(2,1,1)
% barcolor='k';
% 
%    bar(iter,abs(Fobj),0.1,barcolor); 
% 
% axis([0 10 0  max(abs(Fobj))*1.1 ])
% grid on;
% title('Objective function x Iterations') 
% ylabel('Objective function')
% xlabel('Iterations')
%%
% figure(1) % plot CPU time at each iteration
% subplot(2,1,2)
% time = time/(1.3*max(t_ref));
% time_ac = time_ac/(1.3*max(t_ref));
% 
% 
% barcolor='k';
% 
%    bar(iter,time_ac,0.1,barcolor);  
%    hold on;
%    bar(iter,time,0.1,'w'); 
% 
% axis([0 10 0  max(abs(time_ac(is0+1:is)))*1.1 ])
% grid on;
% title('Normalized CPU time x Iterations') 
% ylabel('Normalized CPU time')
% xlabel('Iterations')
% legend('CPU Accumulated','CPU on iteration',1);
%  EMS calculation
k=1;
figure(2) % plot EMS at each iteration  
    barcolor='k';
    for ic=1:nUt 
%          for i=is0+1:is-1
             if plot==1
                MSE(ic,1) = (norm(fit(end,:,ic)-u_opt(end,:,ic)))^2/100;
             else
                MSE(ic,1) = (norm(fitt(end,:,ic)-u_opt(end,:,ic)))^2/100;
             end
                SNR(ic,is) = 10*log10( (norm(u_opt(end,:,ic))^2)/100/MSE(ic,1) );
        %        [PSNR(ic,is),MSE(ic,is),MAXERR(ic,is),L2RAT(ic,is)] = measerr( u_opt(end,:,ic),fit(end,:,ic) );
%          end
             %   out_SPAIDER.psnr = PSNR;
                out_SPAIDER.psnr = 1;
               subplot(2,1,1)
%                bar(iter(1:length(MSE(ic,:))),MSE(ic,:),0.1,barcolor);  
               if nUt==1
               title('MSE x Iterations of variable')     
               else
               title(['MSE x Iterations of control variable',int2str(ic)]) 
               end
               ylabel('MSE')
               xlabel('Iterations')          
%                axis([0 10 0  max(MSE(ic,(is0+1:is)))*1.1 ])
               grid on;
%                subplot(2,1,2)
%                bar(iter,SNR(ic,:),0.1,barcolor);  
               if nUt==1
               title('SNR x Iterations')    
                    else
               title(['SNR x Iterations of control variable',int2str(ic)]) 
               end
               ylabel('SNR')
               xlabel('Iterations')          
%                axis([0 10 0  max(SNR(ic,(is0+1:is)))*1.1 ])
               grid on;
    end


 
%    grid on;
    
    end
end
   if sc==1  && strcmp(optSPAIDER.wavelets.procedure,'on')
        return
   end
 
end


 
  
  

 
 
 