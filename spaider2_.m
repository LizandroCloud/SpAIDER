function[wopt, erropt, out_jac,out_hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv,detais00,...
    detais11,uth,ld0,ld1,md00,md11,vd00,vd11,difU,difd,mseU,is,condh] =  spaider2_(optSPAIDER,optNLP,optODE)
 % ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Química - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 %                              Description
 % This code is applyed to solve dynamic optimization problems with
 % wavelets basis as adaptive startegy for mesh discretization. In 
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
 %*************************************************************************
 
 if strcmp(optSPAIDER.solver.type,'external')
 [ dx y par e_var a_var c_var J cr_eq cr_deq optSPAIDER.input] = dae_model;
 else
 [ dx y par e_var a_var c_var J cr_eq cr_deq optSPAIDER.input] = dae_model;
 end
 [t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,frozen,...
    solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,tunn,fref,typ,wfact,reg,est] = spaider_head(optSPAIDER);
 
 
 w01=[]; % pre-allocating
 w02=[]; detais00=[]; detais11=[]; uth=[]; ld0=[]; ld1=[]; md00=[]; md11=[]; difU=[]; difd=[]; mseU=[]; vd00=[]; vd11=[]; ts=[];
 
 li=0;
 ci=0;
 w0ns=[];
 j =0;
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
 if exist('alg.m') 
 else
  [s]=write_alg('', '', '', 1, 1, 1, 0, ''); % algebraic equations
 end
 if strcmp(optSPAIDER.solver.name,'emso')
     optSPAIDER.solver.type = 'external';
     % do not forget option fr jacobian
 else
%      optSPAIDER.solver.type = 'internal';
 end
  if strcmp(optSPAIDER.solver.type,'external')
%        [s]=write_alg('', '', '', 1, 1, 1, 0, ''); % algebraic equations
    else
 % dynamic functions
 % ************************************************************************
 [param eq eq_sens n_p n_c n_a y J dJ cr_eq d_cr_eq cr_deq d_cr_deq cr_eq_type cr_deq_type] = dyn; % calculating Jacobian
 
 % writing program routines that depens of Jacobian calculation
  [s]=write_alg(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % models
 [s]=write_models(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % models
 [s]=write_obj(J, dJ,param); % objective
 [s]=write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints
 
 %*************************************************************************
  end
 



%  % dynamic functions
%  % ************************************************************************
%  [param eq eq_sens n_p n_c n_a y J dJ cr_eq d_cr_eq cr_deq d_cr_deq cr_eq_type cr_deq_type] = dyn; % calculating Jacobian
%  
%  % writing program routines that depens of Jacobian calculation
%  
%  [s]=write_models(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % models
%  [s]=write_obj(J, dJ,param); % objective
%  [s]=write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints
%  
%  %*************************************************************************
 if optSPAIDER.wavelets.Wav == 1
 if adap == 1
     file = optSPAIDER.wavelets.tptr;
 elseif adap ==0
     file = 'fixed threshold';
 elseif adap==2
     file = 'norm threshold';
 end
 else
     file = 'no-wavelets';
 end
%**************************************************************************

% for mulptible shooting
if strcmp(optSPAIDER.problem.type,'multiple shooting')
    
    
    
    if size(optSPAIDER.input.u_lb)<= size(optSPAIDER.input.x0)
                                                                   
        adv{1}= ' For multiple shooting you must specify the constraints of initial values';
        adv{2}= ' Please, go to routine: model_input.m and specify data.u_lb and data.u_ub';
        status = 'error';
        [mess] = spaider_message(adv,status);
        if mess == 1
            return;
        else
            
        end
%         function [mess] = spaider_message(adv)
%                 ast = repmat('                   *',lengh(adv),1);
%                 disp('                   ******************************************************************************');
%                 disp('                   *                              Warning:                                      *');
%                 disp('                   *                                                                            *');
%                 for kj = 1:2;
%                     message = strcat((ast(kj)),adv{kj'},'   *'); 
%                     disp(message); 
%                 end
%                 disp('                   ******************************************************************************'); 
%         end

    end
    
    
    if typ==1  % if time free
    nUt = nU;  % original number of control variables
    nU = nUt + length(optSPAIDER.input.x0) + 1;  % number of controls + initials + time
    else
    nUt = nU;  % number of controls 
    nU = nUt + length(optSPAIDER.input.x0);  % number of controls + initials 
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
    if strcmp(optSPAIDER.problem.type,'multiple shooting')
        
        
        
        for j = 1+nUt:nUt+nU-1; % for each control, variable
%         ns(j) = (2^is); % initial number of stages
        wL(j,1:1) = u_lb(j-nUt+1)*ones(1,1);  % lower boundary of u
        wU(j,1:1) = u_ub(j-nUt+1)*ones(1,1);  % upper boundary of u
        w0(j,1:1) = u0(j-nUt+1)*ones(1,1);  % initial estimative da u...
%         ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
%         w_frozen = w0;
%         check = [2 3 4 5 6 7 8 9];
%         ts_frozen = ts;
%         ts_a = ts;
        end     
    end
        
        
    end
%  iu=is+2; % graphical counter....
    ns_a = ns(1:nUt);
%     if typ == 1
%             wL(nU,1) =  u_lb(nU);
%             wU(nU,1) =  u_ub(nU);
%             if is==is0
%             w0(nU,1) = u0(nU);  % initial estimative da u...\
%             else
%             w0(nU,1) =  talf;
%             end         
%     end

    disp(sprintf('Iteration: %d', is));  % writing on the terminal
    disp(sprintf('Number of stages: %d', ns));
    disp(sprintf('Method: %s',file)) 
    % calling optimization algorithm...
     % pre-optimization#####################################################
    k=0;
     tss(1) = ts(1,1);
	ki=1;
    w00=ones(1,nU*1);
    tss=[];
    for i=1:nUt
         for j=k+1:ns_a(i)+k
             tss(j)=ts(i,ki);
             ki=ki+1;
         end
         tss(j+1)=ts(i,ki);
         k = j+1;
         ki=1;
    end
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
 ko=k;
 ts_a = ts_frozen;
 
 if strcmp(optSPAIDER.problem.type,'multiple shooting')
        k=0;
        j=1:nU-1
        w00(j+ko) = w0(j+1,1);
        wLL(j+ko) = wL(j+1,1);
        wUU(j+ko) = wU(j+1,1);
        
 end
 
 if typ==1
 nss = sum(ns)-typ;
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
 


 % optimization
 %########################################################################
%  [f_inicial] = objfunc(x0,nss,tss,w00a,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER)

% calling optimization routine...
 tic
 ob = @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER);
 [wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER), w00, [], [], [],...
 [], wLL, wUU,  @(ws)constr(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER), optNLP);  
%  [wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER), w00, [], [], [],...
%  [], wLL, wUU,  [], optNLP);  
 toc
 
 % storaging hessian...
 szh=size(hessian);
 condh(is) = cond(hessian); 
 szg=length(grad);
 out_hessian(is,1:szh(1),1:szh(1)) = hessian;
 out_jac(is,1:szg) = grad;
 
 if is==9
     klj=is;
 end
 
 iter(is)=is;
 time(is)=toc;
 time_ac(is)=time(is)+toc_ant;
 toc_ant=time_ac(is);
 kit=kit+1;
 
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


% [f_final] = objfunc(x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER)
 
 % pos optimization********************************************************
 
 if typ==1
    talf=wopt(end);
 else
    talf=1;
 end

% State graphs
if is>is0
[r0 tst fst] = g_state(plot_state,nss,tss,wopt,ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,talf,reg,est,optSPAIDER);
end        

%############ reconfigurating wopt for frozen option 
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
    wopt=[];
    wopt=wxi;
    ts_a=ts;
    end;
%#########################################################
    ic = 1:nUt;
    wx(ic,1:length(wopt(ic,:)),is)=wopt(ic,:);  % each line of w is a control variable...
    tx(ic,1:length(ts(ic,:)),is)=ts(ic,:); % time stages
    nsx(ic,1,is)=ns(ic);  % number of stages

for ic = 1:nUt  % nUt is the number of control variables
%      if ic==1
%          tptr = 'sqtwolog'; scal = 'sln'; % scale procedure
%      end
%      if ic==3
%          tptr = 'sqtwolog'; scal = 'sln'; % scale procedure
%      end
if is==8
    aq=1;
end
    if opC==1 % Wavelets analisys
        if (is>is0) && (is<=isF)    % if the iteration is superior than first iteration and less then last iteration 
            wopt(ic,1:(ns(ic))) = wopt(ic,1:(ns(ic)))'; %transposing
            [cfd,detais0,maxlev,C,L] = wav_mesh(wopt(ic,1:(ns(ic))),is,ic,wname); % into wavelets domain
            detais00(ic,1:length(C),is-is0)=C;  % details before thresholding
            ld0(ic,is-is0)=length(C);
            [md0,vd0,muci,sci] = normfit(C);
            md00(ic,is-is0)=md0;
            vd00(ic,is-is0)=vd0;
            [Cth, thrs1, detais1, difdd, w01, s, threshold] = wav_thresh(C,L,maxlev,wname,detais0,wopt,tptr,sorh,adap,scal,ic,is,is0,optSPAIDER.wavelets.normfrac,optSPAIDER); % Wavelets thresholding
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
            [Pgrad] = adapt(w01,Cth,L,maxlev,ic,tunn,ts,adap,optSPAIDER.wavelets.normfrac) % prospective points...
            
            
            Pgrad=Pgrad(Pgrad~=0);
%             if output.iterations < 10 
%                 Pgrad = downsample(Pgrad,2)
%             end
%           stop barrier
            PgradX(ic)=norm(Pgrad);
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
                    return
                elseif   stop_flag_norm ==0
                    avc=1
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
                [r_ix2a] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is);
                % Control profile before thresholding
                [r_ix3a] = g_ix3a(ix3,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is);
                % Details bars before thresholding
                [r_ix2b] = g_ix2b(ix2,ic,ns_a,talf,C,nUt,L);
                 % Details bars after thresholding
                [r_ix3b] = g_ix3b(ix3,ic,ns_a,talf,C,nUt,L);
                % Thresholded control
                [r_ix4a fit] = g_ix4a(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,Cth,nUt,w01,wopt,is); 
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
                [r_ix2a] = g_ix2a(ix2,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is);
                % Control profile before thresholding
                [r_ix3a] = g_ix3a(ix3,ic,ns,ns_a,t0,tf_p,ts,talf,C,nUt,wopt,is);
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
% grid refinement  (meshref function)...
%##########################################################################

[ts_ref w_ref ns_ref point] = meshref(Pgrad,ns,ic,ts_a,w02,t0,tf,wfact);
ts(ic,1:ns_ref+1)=ts_ref;
w0ns(ic,1:ns_ref)=w_ref;
ns(ic)=ns_ref;  
 
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
end
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
 w0ns=[];
 wLL=[];
 wUU=[];
 ts_r =[];
 
 % Stop Criterium 
 if is>is0
 if opC==1

 
 [sc cvc cdc] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,u_opt,cvc1,cdc1,time_ac,nUt);
  cvc1=cvc;
  cdc1=cdc;
   conv(is)=cvc1(end);
   convv(is)=cdc1(end);
   if sc==1
       return
   end
 
 else
   [sc cvc cdc] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,u_opt,cvc1,cdc1,time_ac,nUt);
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
%%  EMS calculation
% k=1;
% figure(2) % plot EMS at each iteration  
%     barcolor='k';
%     for ic=1:nUt 
         for i=is0+1:is-1
             if plot==1
                MSE(ic,i) = (norm(fit(i,:,ic)-u_opt(i,:,ic)))^2/100;
             else
                MSE(ic,i) = (norm(fitt(i,:,ic)-u_opt(i,:,ic)))^2/100;
             end
%               SNR(ic,i) = 10*log10( (norm(u_opt(i,:,ic))^2)/100/MSE(ic,i) );
         end
%                subplot(2,1,1)
%                bar(iter(1:length(MSE(ic,:))),MSE(ic,:),0.1,barcolor);  
%                if nUt==1
%                title('MSE x Iterations of variable')     
%                else
%                title(['MSE x Iterations of control variable',int2str(ic)]) 
%                end
%                ylabel('MSE')
%                xlabel('Iterations')          
%                axis([0 10 0  max(MSE(ic,(is0+1:is)))*1.1 ])
%                grid on;
%                subplot(2,1,2)
%                bar(iter,SNR(ic,:),0.1,barcolor);  
%                if nUt==1
%                title('SNR x Iterations')    
%                     else
%                title(['SNR x Iterations of control variable',int2str(ic)]) 
%                end
%                ylabel('SNR')
%                xlabel('Iterations')          
%                axis([0 10 0  max(SNR(ic,(is0+1:is)))*1.1 ])
%                grid on;
%     end


 
%    grid on;
    
    end
end
end


 
  
  

 
 
 