function [ f, t_ult, t_out, y, ti, df ] = function_ode( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER )
% global reg;
% global est;
% global DEBUG;  % for dassl solver...
% global t_ult;
global modelo_dFdx modelo_dFdu
%% Parameters
% DEBUG  = 1; 

dae_index=0;
t0 = 0;
ks=1;
d_ode = optSPAIDER.solver.sampling;
% d_ode = ne(1)/2;
% nearestEven = floor(d_ode) - mod(floor(d_ode), 2);
% d_ode = nearestEven;
yi = [];
% reg = 0; % regularization...
% fix =0;  % fix solver...
est = 1;  % integration by steps...
ti = zeros(1,500);
ts = roundn (ts, -5);
koz=1;  %reorganizing time steps...
k=1;
ko=0;
ne_a=ne;
rtol= 1e-10;       % relative tolerance for dassl
atol= 1e-10;       % absolute tolerance for dassl



dim=size(ws);
ts_ms = round(ts);
ne_ms=(dim(2)-n_c*sum(ne))/n_s -0; % numer of elementos of multiple shooting
if strcmp(optSPAIDER.problem.method,'multiple shooting')  % keeping the initial conditions for multiple shooting
    for j=1:n_s
    x0_ms(1:ne_ms,j) = ws(sum(ne) + (j-1)*ne_ms + 1 : sum(ne) + (j)*ne_ms);
    end
end

x(3) = ws(1);
 if typ==1
%       tf=ws(end);
%     if isempty(t_ult)
%         t_ult=tf;
%     end
   tal = ws(end);
%     for i=1:length(ts)
%         ts(i)=ts(i)*ws(end)/t_ult;
%     end
 else
     tal=1;
 end

%% Organizing vectors 
kq=1;
for ix=1:n_c
    for j=1:ne(ix)
        w0(ix,j)=ws(kq);
        kq=kq+1;
    end
end

if frozen==1
    ne=ns_a;
end

for ix=1:n_c
    for j=1:ne(ix)+1
    te(ix,j)=ts(koz);
    koz=koz+1;
    end
end

if frozen==0                
         for ix=1:n_c
        for j=1:ne(ix)
%         for k = 1:length(check(i,:))
%             if j == check(i,j)
%                   wx(i,j)=ws(ko);
%                   ko=ko+1;
%             else
                  wx(ix,j)=ws(ko+1);
                  ko=ko+1;
%             end
        end
        
    end
else

    if is == is0
     for ix=1:n_c
        for j=1:ne(ix)
%         for k = 1:length(check(i,:))
%             if j == check(i,j)
%                   wx(i,j)=ws(ko);
%                   ko=ko+1;
%             else
                  wx(ix,j)=ws(ko+1);
                  ko=ko+1;
%             end
        end
        
    end
else
    for ix=1:n_c
    for j=1:ne(ix)
        for k = 1:length(check(ix,:))
            if j == check(ix,k)
                  wx(ix,j)=w0(ix,k);
                  
                  break;
                 
            else
                  wx(ix,j)=w_froz(ix,j);
            end
            
        end
  
    end
    end
    end
    end
% ne=ne_a;
ts_frozen = roundn (ts_frozen, -5);
te = roundn (te, -5);
t_out = [];
zssx= []; 


%% Forward state integration with numerical jacobian...
if nargout ~= 6 % Forward state integration with numerical jacobian...

zss= [];   

z0 = [ x0 ];
 if frozen==0
    check = 2:length(ts);
    for ix=1:n_c
        for j=2:ne(ix)
            check(ix,j) = j;
        end
    end
 end 
 %=========================================================================
 
 
     ko=0;
     if est==1  % with integration by steps or elements...
        
        if strcmp(optSPAIDER.problem.method,'multiple shooting') 
                           ti = unique(abs(ts));
%                            ti=ti/60;
                 %         ti =num2cell(ti);
                 %         te =num2cell(te);
                 tix =  length(ti)-1;
                 wxx=length(wx(:,1))   ; 
%         ti = unique(abs(ts));  %organizing time steps...
                ug=  num2cell(wx);
%         ne=num2cell(ne);
                disc=10;
%                 tic
                for (ks = 1:tix)  % parallel
                % for (ks = 1:tix)  % serial
                        u={}; nex={}; ugx={}; tex={};
                        [to tf dtt] = f_dt(ti,ks,disc);
%                 tex{ks,:}=te{1,:};
%                 ai =ti{ks};
%                 z0 =  x0(2:end) ;
%                  i=1;
%                 ti = unique(abs(ts));  %organizing time steps...
%                 u = zeros(wxx);
%                         pzo = {z0};
                        pzo = {f_zo(ks,{z0}, {x0_ms})};
                        for i = 1:wxx 
                            nex = nstages(ne,ks,i);
                            ugx = f_ug(i,wx);
                            [tex tex_o tex_1] = f_te(i,te,ks);
                            
%                       u{i} = u_reg1(0.5*(ti{2}+ ti{3}),[te{i,:}],[ug{i,:}],[(ks)],[ne{i}]);  % calculating the value of u in time domain (ti)
                            u{ks,i} = u_reg1(0.5*(tex_o+ tex_1),tex,ugx,1,nex);
%                      u{ks,i}=ug{i};
                        end

%                         [tspan,hh,zs]= parsim('emso',[0 :dtt: tf-to],[  tf-to]', [u{ks,:} pzo]); % original
                    [tspan,zs]  = optSPAIDER.solver.func(@(t,x)modelo_eq(t,x,cell2mat(u)',(ks),ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
                        tal,est,optSPAIDER.prob.uref,...
                     optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], cell2mat(pzo), optODE );
                    % h = 30;
                    % h=(ti(ks+1)-ti(ks))/(ns*10);
                    % [tspan,zs]  = quadrature_collocation(@(t,x)modelo_eq(t,x), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], cell2mat(pzo), h, cell2mat(u)', ne,te,n_s,n_c,ne_a,ts_frozen,check,...
                         % tal,est,optSPAIDER.prob.uref,...
                      % optSPAIDER.prob.tref,optSPAIDER.prob.nref);
%  
%                          [tspan,hh,zs]= parsim('emso',[0 :0.5: 1],[  tf-to]', [u{ks,:} dtt pzo]); % original
%                           [tspan,hh,zs] = sim('emso',[ti(1) :0.5: 1],[],[(ti(ks+1)-ti(ks))' u (ti(ks+1)-ti(ks)) ones(1,1)*z0 ]);
%                        [tspan,hh,zs] = sim('emso',[ti(1) :0.5: 1],[],[(ti(ks+1)-ti(ks))' u (ti(ks+1)-ti(ks)) ones(1,1)*z0 ]);
%                [tspan,hh,zs] = sim('emso',[ti(ks) :( (ti(ks+1)-ti(ks))/10 ): ti(ks+1)],[],[ti(ks+1)' u ones(1,1)*z0]);
                 %[tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks) ti(ks+1)], z0, optODE );
                         tspan = tspan*1;
                        tsc{ks} = tspan;
                        zsc{ks} = zs;
%                  for ivf = 1+ko:length(tspan)+ko
%                      [yout] = alg_eq(zs(ivf-ko,:),u);
%                      y(ivf,:) = yout;
%                  end
%                  ko = length(tspan)+ko;
%                  t_out= [ t_out; tspan];
%                  z0 = zs(end,:);
%                  zss= [zss;zs];
                end
                                 zss=[]; t_out=[];
     for kx=1:tix
         zsa =  (zsc{1,kx});
         zss = [zss;zsa];
         tspan = (tsc{kx});
%          t_out= [ t_out; (kx-1)*tspan(end)+tspan];
         t_out= [ t_out; tspan(1:end)'];
     end

     for ivf = 1:tix
                for i = 1:length(wx(:,1))  % length of each control vector...
                  ui(ivf,i) = u_reg1(1*(ti(ivf)+ ti(ivf+1)),te(i,:),wx(i,:),1,ne(i)); % u_reg1 is a function heaviside regularized...
                end
         [yout] = alg_eq(zss(ivf,:),ui);
         y(ivf,:) = yout;
     end

     else  % for single shooting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      A='LIT_004.SPinput -'; 
%                      B=num2str(wx(1,1));
%                      C=[A,B];
%                      [pfd, y0,ydot0] = UOTR(C);
%                      z0 = y0; 
        ti = unique(abs(ts));  %organizing time steps...
%         z0=y0;
        tic
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
%                 ts(end)
                    
%                   pAD = myAD(u);
%                  options = CVodeSetOptions('RelTol',1.e-8,...
%                           'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
%                           'LinearSolver','Dense',...
%                           'JacobianFn',@djacfn,...
%                           'RootsFn',@rootfn, 'NumRoots',2);
                      
%                 CVodeInit(@rhsfn, 'BDF', 'Newton',ti(ks), z0, options);
                
%                t=0;
                     
%                  [zs]  = ode5(@modelo_eq, [ti(ks):0.001: ti(ks+1)], z0,u' );
%                  global UX
%                  UX=u(i);
%                  z0(210)=u(i);
                 
                 
%                     [tspan,zs] = ode15i('daesys',[ti(ks):(ti(ks+1)-ti(ks))/50:ti(ks+1)],z0,ydot0);
%                  if length(tspan)<10
%                      warning('Falha na Integra��o. Gerando novas condi��es iniciais consistentes');
%                      A='LIT_004.SPinput -'; 
%                      B=num2str(UX);
%                      C=[A,B];
%                      [pfd, y0b,ydot0] = UOTR(C);
%                      [tspan,zs] = ode15i('daesys',[ti(ks):(ti(ks+1)-ti(ks))/50:ti(ks+1)],z0,ydot0);
%                  end
                    %ORIGINAL
                    for j=1:n_c
                        yi{j}=1;
                    end
                
                 [tspan,zs]  = optSPAIDER.solver.func(@(t,x)modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,yi,optSPAIDER.prob.uref,...
                     optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], z0, optODE );
%                     [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
%uws,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks):(ti(ks+1)-ti(ks))/50:ti(ks+1)], z0, optODE );
%                 ydot0(210)=0;
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 
                 if ks==1
                    ko = length(tspan)+ko;
                    t_out= [ t_out; tspan];
                    z0 = zs(end,:);
                    zss= [zss;zs];
                 else
                    ko = length(tspan)+ko;
                    t_out= [ t_out; tspan(2:end)];
                    z0 = zs(end,:);
                    zss= [zss;zs(2:end,:)]; 
                 end
                 
                 
%                  ko = length(tspan)+ko;
%                  t_out= [ t_out; tspan];
%                  z0 = zs(end,:);
%                  zss= [zss;zs];
        end
        % disp(toc)
        end
     else   % if est=0, full integration
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(ks,i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
        end
%            [tspan,zs]  = ode23( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
            for j=1:n_c
                % for i = 1:ne_a(j) % a priori interpolation if est=0...
                    x = te(j,2:ne_a(j)+1); 
                    yx = u(1:ne_a(j),j)';
                    yi{j} = interp1(x,yx,optSPAIDER.problem.interpolation, 'pp');  % interpolacao ativada
%                  u(i)= ppval(yi,t);
                % end
            end
            span = [0:ts(end)/d_ode: ts(end)];
            tic
            [tspan,zs]  = optSPAIDER.solver.func( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,...
            check,tal,est,yi,optSPAIDER.prob.uref,optSPAIDER.prob.tref,ws),...
            span,z0, optODE );
            % disp(toc);
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 
 
                 
                 ko = length(tspan)+ko;
            t_out= [ t_out; tspan];
            z0 = zs(end,:);
            zss= [zss;zs];

            % s = modelo_eqsens_nonlin(t_out,x0,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,[],ws,tal,est,yi,[],[],[],modelo_dFdx,modelo_dFdu);
     end
%      elseif    strcmp(solver,'ode23s')   
%      if est==1  % with integration by steps...
%         ti = unique(abs(ts));  %organizing time steps...
%         for ks = 1:length(ti)-1
%                 for i = 1:length(wx(:,1)) 
%                   u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
%                 end
% %                 ts(end)
%                  [tspan,zs]  = ode23s( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est), [ti(ks) ti(ks+1)], z0, optODE );
%                  t_out= [ t_out; tspan];
%                  z0 = zs(end,:);
%                  zss= [zss;zs];
%         end
%      else   % if est=0, full integration
% %            [tspan,zs]  = ode23( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
%   
%            [tspan,zs]  = ode23s( @(t,x) modelo_eq(t,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est), [ts(ks) ts(ns+1)], z0, optODE );
%             t_out= [ t_out; tspan];
%             z0 = zs(end,:);
%             zss= [zss;zs];
%      end
 
 %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 f = zss(:,:);
 if est==1
 % atime=linspace(0,t_out(end),d_ode*ns_a + 1);
 else
 % atime=linspace(0,t_out(end),d_ode*1 + 1); 
 end
 % f =  interp1(t_out,zss,atime,'linear');
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
 
%%
else
%% Forward state integration with Analytic (user specified) Jacobian...
  z0_emso = x0; % only for emso simulator
  zss= [];
  if strcmp(optSPAIDER.problem.method,'multiple shooting') 
     zdf = ones(n_s*(ns)+n_s*(ns-1),1);
  else
     zdf = ones(n_s*(ns),1);
  end
  zdf = 0*zdf;
  % z0 = [ x0' ; zdf ];
  z0 = [ x0' ; zdf ];
  z0x = z0;
  if frozen==0
%     check = 2:length(ts);
    for i=1:n_c
        for j=2:ne(i)
            check(i,j) = j;
        end
    end
 end 
 
 
   %#####################################################################
 % ODE45 MODEL
 %#####################################################################
 Rz0 = [ x0 ];
 if reg==0
     ti = unique(abs(ts));  %organizing time steps...
     ko=0;
            for ks = 1:length(ti)-1
%              ko=0;   
             
                for i = 1:n_c 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end

                %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                for j=1:n_c  % analysing if time corresponds to the actual control stage...
                     for i = 1:ne(j) 
                        tx = [0.5*(ti(ks)+ ti(ks+1)) ts_frozen(j,2:ne(j))];
                        tx = sort(tx);
                        posi = find(tx == 0.5*(ti(ks)+ ti(ks+1)),1,'last');
                
                
                        if ((posi)== i)
                            uws(j,i) = 1;
                            else   % in contrary...
                            uws(j,i) = 0;
                        end
                     end
                end
                 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
%             for i = 1:length(ws(:,1)) % a priori interpolation if est=0...
%                 x = ts(i,:); 
%                 yx = ws(i,:);
%                 yi = interp1(x,[yx yx(end)],optSPAIDER.problem.interpolation, 'pp');  % interpolacao ativada
% %                 u(i)= ppval(yi,t);
%             end

                 % pseudo_sen
                 
                 [tspan,zs]  = optSPAIDER.solver.func(@(t,x)modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,yi,optSPAIDER.prob.uref,...
                     optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], z0(1:n_s), optODE );


                % Sensibilidade Analitica...

                % [~,zsx]  = optSPAIDER.solver.func( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,uws,tal,...
                    % est,yi, optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], z0x, optODE );
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 
                 if ks==1
                    ko = length(tspan)+ko;
                    t_out= [ t_out; tspan];
                    z0 = zs(end,:);
                    % z0x = zsx(end,:);
                    zss= [zss;zs];
                    % zssx= [zssx;zsx];
                 else
                    ko = length(tspan)+ko;
                    t_out= [ t_out; tspan(2:end)];
                    z0 = zs(end,:);
                    % z0x = zsx(end,:);
                    zss= [zss;zs(2:end,:)];
                    % zssx= [zssx;zsx(2:end,:)];
                 end
                 zxp = zss;
                 
%                 ko = length(tspan)+ko;
%                 t_out= [ t_out; tspan];
%                 z0 = zs(end,:)';                
%                 zss= [zss;zs];
%                  zss = interp1(t_out,zss,0:.1:6);
            end



%   zss = interp1(t_out,zss,0:.1:6);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Sensitivity analytical solution with variable jacobian
    % if flagstart    %No need to update jacobians because ws(1:Nts)=ws(1)
        % dFdx(:,:,1)=modelo_dFdx;
        % dFdu(:,:,1)=modelo_dFdu;
    
        % s=modelo_eqsens_ctrl_Kp(t_out,s(1,:),ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx(:,:,1),dFdu(:,:,1));
    
        % flagstart=false;
    % else
        sv=[];
        for j=1:n_c
            Id=zeros(n_s,ns_a(j));
            sv=[sv';Id(:)]';
        end
        s(1,:)=sv;
        Nts = length(ti)-1;
        
        idts=zeros(Nts,1);
        for k=1:Nts
            idts(k)=max(find(ts(k)>=t_out));
        end

        % wz = rpmat(ws,)

        for k=1:Nts


                for i = 1:n_c 
                  wsu(i) = u_reg1(0.5*(ti(k)+ ti(k+1)),te(i,:),wx(i,:),k,ne(i)); 
                end
                % wsu=wsu';
            % if k>1
                %Update dFdx
                fhand=@(t,x)modelo_eq(t,x,wsu',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,yi,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref);

                % fhand=@(t,x)modelo_eq(t,x,wsu(:,k),ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,yi,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref);
                            % modelo_eq(t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
                fval=fhand([],zss(k,:)');
                [dFdx(:,:,k),~,~]=numjac(fhand,[],zss(k,:)',fval,zeros(n_s,1),[],0);
                dFdx(isnan(dFdx)) = 1e-5;
                %Update dFdu
                wsu(:,find(wsu==0))=eps;
                fhand=@(t,u)modelo_eq(t,zss(k,:)',u,ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,yi,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref);
                % fhand=@(t,u)uotr_ctrl_Kp(t,zss(ti(k),:)',[u(1);uf(k,2:15)';u(2:length(idu))],1,[],modelo_par);
                % fval=fhand([],wsu(:,k));
                fval=fhand([],wsu');
                % [dFdu(:,:,k),~,~]=numjac(fhand,[],wsu(:,k),fval,zeros(n_c,1),[],0);
                [dFdu(:,:,k),~,~]=numjac(fhand,[],wsu',fval,zeros(n_c,1),[],0);
            % else
                % dFdx(:,:,k)=modelo_dFdx;
                % dFdu(:,:,k)=modelo_dFdu;
            % end        
            
                % if k<Nts %From ts(k) to ts(k+1)
                % sv=p_sens( [ti(k):(ti(k+1)-ti(k))/d_ode: ti(k+1)] ,s(end,:),ws,ks,ti,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,check,tal,est,k,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref,dFdx(:,:,k),dFdu(:,:,k));
                sv=p_sens([ti(k):(ti(k+1)-ti(k))/(d_ode): ti(k+1)],  s(end,:),ks,te,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,check,tal,est,k,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref,dFdx(:,:,k),dFdu(:,:,k));
                % s=[s(1:end,:);sv(2:end,:)];
                s=[s(1:end,:);sv(2:end,:)];
                % else %From ts(end) to t_out(end)
                % if idts(end)<length(t_out)
                    % sv=modelo_eqsens_ctrl_Kp((t_out(idts(k)):1:t_out(end))',s(1,:),ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx(:,:,k),dFdu(:,:,k));
                    % sv=p_sens([ti(k):(ti(k+1)-ti(k))/(d_ode): ti(end)],s(end,:),ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,check,tal,est,k,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref,dFdx(:,:,k),dFdu(:,:,k));
                    % s=[s(1:end-1,:);sv((1:tstep:(t_out(end)-t_out(idts(k))+1))',:)];
                    % s=[s(1:end,:);sv(2:end,:)];
                % end
                % end
        end
    % end
    
    %Stocking sensitivities into block vector in f_out
    zss=[zss,s];
    % disp(zss);
    % zss(:,1:3) = zxp(:,1:3);
    % zss = zssx;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Axx=zss-zssx;
    % disp( (Axx) );


  %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 if strcmp(optSPAIDER.problem.method,'multiple shooting')
 for is = 1:(ns)
    df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 ki=is; nimax=(2*ns - 1);
  for is = ki:ki
    df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 
 
 else

 for is = 1:ns
    df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 end
 f = zss(:,:);
 % pass = t_out(end)/d_ode;
 % atime=linspace(0,t_out(end),d_ode*ns_a + 1);
 % f =  interp1(t_out,zss,atime, 'spline');
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
%      optSPAIDER.problem.analitic_derivarive='automatic differentiation';
     if strcmp(optSPAIDER.problem.analitic_derivarive,'automatic differentiation' )
     t_out=[]; zsa = [(zeros(3,n_s*(length(ti)-1)))]; kt=0;
      z0 = [x0'; (zeros(1,n_s))'];
      ko=0;
     if est==1  % with integration by steps or elements...
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
%                 ts(end)
%                   pAD = myAD(u,x);
                 [tspan,zs]  = optSPAIDER.solver.func(@(t,x)modelo_eq(t,x,u' ,ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,...
                     optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/d_ode: ti(ks+1)], z0, optODE );
%                    [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
%uws,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks):(ti(ks+1)-ti(ks))/50:ti(ks+1)], z0, optODE );

                     ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 
                      kv = 1:n_s
                         zsa(kt+1:length( tspan)+kt, 1:n_s) = zs(:, 1:n_s);
                         zsa(kt+1:length( tspan)+kt, (n_s*(ks)+1):n_s*(ks+1)) = zs(:, n_s+1:n_s*2);
%                      end
                 
                 
                 ko = length(tspan)+ko;
                 t_out= [ t_out; tspan];
                 z0 = zs(end,:);
                 kt = length( tspan)+kt;
%                  zss= [zss ;zsa];
%                  zsa = [(zeros(3,n_s*(length(ti)-1)))];
            zssc = zsa;
        end
        
        
        
     else   % if est=0, full integration
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(ks,i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
        end
%            [tspan,zs]  = ode23( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
            [tspan,zs]  = optSPAIDER.solver.func( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [0 ts(ns+1)], z0, optODE );
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
            t_out= [ t_out; tspan];
            z0 = zs(end,:);
            zss= [zss;zs];
%               zss = interp1(t_out,zss,0:.1:6);
     end
 
     
   %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 for is = 1:ns
 dfc(:,1:n_s,is) = zssc(:,is*n_s+1:is*n_s+n_s);
 end
 f = zss(:,:);
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
     end
end
 

end









