function[wopt, erropt, grad,hessian,wx,tx,nsx,time,time_ac,Fobj,thrs,conv,convv,detais0,detais1] =  spaider3_(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,opC,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,optNLP,optODE,tunn,fref,typ,wfact,reg,est,optCONV,optWAV,optSPAIDER,optPLOT)
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
 %*************************************************************************
 
 w01=[]; % pre-allocating
 w02=[];
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
 iter=[1 2 3 4 5 6 7 8 9 10];
 Fobj=zeros(1,10);
 time=zeros(1,10);
 time_ac=zeros(1,10);
 MSE=zeros(1,10);
 SNR=zeros(1,10);
 toc_ant=0; 
 conv(1) =10;
 cvc1=[0];
 cdc1=[0];
 if strcmp(optSPAIDER.solver.name,'emso')
     
     optSPAIDER.solver.type = 'external';
 else
     optSPAIDER.solver.type = 'internal';
 end
 if strcmp(optSPAIDER.solver.type,'external')
 else
 % dynamic functions
 % ************************************************************************
 [param eq eq_sens n_p n_c n_a y J dJ cr_eq d_cr_eq cr_deq d_cr_deq cr_eq_type cr_deq_type] = dyn; % calculating Jacobian
 
 % writing program routines that depens of Jacobian calculation
 
 [s]=write_models(eq, eq_sens, param, n_p, length(x0), n_c, n_a, y); % models
 [s]=write_obj(J, dJ,param); % objective
 [s]=write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type); % constraints
 
 %*************************************************************************
 end
 if adap == 1
     file = 'wavelets threshold';
 else
     file = 'heuristic threshold';
 end
%**************************************************************************
 
if typ == 1   % flag that specifies if is time free or not
    nUt = nU-1;
else
    nUt=nU;
end 
 
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
         w00(j)=(w0(i,ki)-wL(i,ki))/(wU(i,ki)-wL(i,ki));
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
 ts_a = ts_frozen;
 if typ==1
 nss = sum(ns)-typ;
 w00(nss+1) = w0(nU,1)-u_lb(nU)/(u_ub(nU)-u_lb(nU));
 wLL(j+1)=0;
 wUU(j+1)=1;
 else
 nss = sum(ns);
 end
 

 
 % optimization
 %########################################################################
%  [f_inicial] = objfunc(x0,nss,tss,w00,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU)

% calling optimization routine...
 tic
[wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER), w00, [], [], [],...
[], wLLn, wUUn, @(ws)constr(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER), optNLP);  
%  [wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)objfunc(x0,nss,tss,ws,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU), w00, [], [], [],...
%  [], wLLn, wUUn,[], optNLP);
 toc

 
 
 iter(is)=is;
 time(is)=toc;
 time_ac(is)=time(is)+toc_ant;
 toc_ant=time_ac(is);
 kit=kit+1;
 
  % correcting normalization
  k=0; ki=1;
  for i=1:nUt
     for j=k+1:ns(i)+k
         wopt(j)=wopt(j)*(wUU(j)-wLL(j))+wLL(j);
         ki=ki+1;
     end
     k = j;
     ki=1;
  end
 
 
%  [f_final] = objfunc(x0,nss,tss,wopt,ns,length(x0),nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est)
 
 % pos optimization********************************************************
 
 if typ==1
    talf=wopt(end)
 else
    talf=1;
 end

% State graphs
r0 = g_state(plot_state,nss,tss,wopt,ns,x0,nUt,w_frozen,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,f,t0,tf,tf_p,talf,reg,est);
        

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
    if opC==1 % Wavelets analisys
        if (is>is0) && (is<isF)    % if the iteration is superior than first iteration and less then last iteration 
            wopt(ic,1:length(ns(ic))) = wopt(ic,1:length(ns(ic)))'; %transposing
            [cfd,detais0,maxlev,C,L] = wav_mesh(wopt,is,ic,wname); % into wavelets domain
            [Cth, thrs1, detais1, difd, w01, s, threshold] = wav_thresh(C,L,maxlev,wname,detais0,wopt,tptr,sorh,adap,scal,ic,is,is0); % Wavelets thresholding
            thrs(ic,1:length(thrs1),is-is0)=thrs1; 
            [Pgrad] = adapt(w01,Cth,L,maxlev,ic,tunn,ts,adap) % prospective points...
                if ( isempty(Pgrad)) 
                return
            end 
            w02(ic,: ) = wopt(ic,: ); % backing to original optimal profile...
            if plot ==1
                ix = (is-is0)*10+ic-1;  % index for pictures generation
                ix2 = ix*10;
                ix3 = ix*100+1;
                ix4 = ix*1000+2;
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
                [r_ix4a fitt] = g_ix4a(ix4,ic,ns,ns_a,t0,tf_p,ts,talf,Cth,nUt,w01,wopt,is); 
                % into wavelets domain
                [cfdth,detais0,maxlev] = wav_mesh(w01,is,ic,wname);
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
        if is>is0  % for first iteration...
            for i = 1:length(wopt(ic,: ))
                wopt(ic,i) = wopt(ic,i) + rand(1)*1e-4;  % forcing swmall variability...
                w02(ic,:)=wopt(ic,: );
            end
        else
            w02(ic,1:ns(ic))=wopt(ic,1:ns(ic));
            w0ns(ic,1:ns(ic) )=wopt(ic,1:ns(ic));
        end
        Pgrad=2:length(wopt(ic,: ))
                convv = 0;
        detais0 = 0;
        detais1 = 0;
        if (is>is0) && (is<=isF) 
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
        % Threshding in bars
        [r_ix3d] = g_ix3d(ix3,is,is0,ic,thrs,L,C);
    else
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
 [sc cvc cdc] = spaider_conv(is,is0,ic,Fobj,f_tol,f_tol2,fref,plot,iter,time,t_ref,fitt,u_opt,cvc1,cdc1,time_ac,nUt);
  cvc1=cvc;
  cdc1=cdc;
   conv(is)=cvc1(end);
   convv(is)=cdc1(end);
   if sc==1
       return
   end
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
%         for i=is0+1:isF-1
%                MSE(ic,i) = (norm(fit(i,:,ic)-u_opt(i,:,ic)))^2/100;
%                SNR(ic,i) = 10*log10( (norm(u_opt(i,:,ic))^2)/100/MSE(ic,i) );
%         end
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


 
  
  

 
 
 