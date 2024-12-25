function [ f, t_ult, t_out, y, ti, df ] = function_dasslc( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER )
 
% global reg;
% global est;
 global DEBUG;  % for dassl solver...
% global t_ult;
global u
% global ts
global ks

%% Parameters
DEBUG  = 1; 
dae_index=0; index = zeros(1,n_s);

t0 = 0;
ks=1;
% reg = 0; % regularization...
% fix =0;  % fix solver...
% est = 1;  % integration by steps...
ti = zeros(1,500);
ts = roundn (ts, -5);
koz=1;  %reorganizing time steps...
k=1;
ko=0;
ne_a=ne;
rtol= 1e-4;       % relative tolerance for dassl
atol= 1e-4;       % absolute tolerance for dassl

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
for i=1:n_c
    for j=1:ne(i)
        w0(i,j)=ws(kq);
        kq=kq+1;
    end
end

if frozen==1
    ne=ns_a;
end

for i=1:n_c
    for j=1:ne(i)+1
    te(i,j)=ts(koz);
    koz=koz+1;
    end
end

if frozen==0                
         for i=1:n_c
        for j=1:ne(i)
%         for k = 1:length(check(i,:))
%             if j == check(i,j)
%                   wx(i,j)=ws(ko);
%                   ko=ko+1;
%             else
                  wx(i,j)=ws(ko+1);
                  ko=ko+1;
%             end
        end
        
    end
else

    if is == is0
     for i=1:n_c
        for j=1:ne(i)
%         for k = 1:length(check(i,:))
%             if j == check(i,j)
%                   wx(i,j)=ws(ko);
%                   ko=ko+1;
%             else
                  wx(i,j)=ws(ko+1);
                  ko=ko+1;
%             end
        end
        
    end
else
    for i=1:n_c
    for j=1:ne(i)
        for k = 1:length(check(i,:))
            if j == check(i,k)
                  wx(i,j)=w0(i,k);
                  
                  break;
                 
            else
                  wx(i,j)=w_froz(i,j);
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
%% Forward state integration with numerical jacobian...
if nargout ~= 6 % Forward state integration with numerical jacobian...

 zss= [];   
z0 = [ x0 ];
 if frozen==0
    check = 2:length(ts);
    for i=1:n_c
        for j=2:ne(i)
            check(i,j) = j;
        end
    end
 end 
 %=========================================================================
 
 
     ko=0;
     if est==1  % with integration by steps or elements...
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
                
                    % calling dasslc...
                    %---------------------------------------------------
                    % [tspan,zs]  = dasslc( @(t,x)state(t,x,ws,ks,ts,reg), [ts(ks) ts(ks+1)], z0, optODE );
                    rpar(1:length(ws),1) = ws(1,:);  % it is a vector....
                    rpar(length(ws)+1,1) = ks;
                    rpar(length(ws)+2:length(ws)+length(ts)+1,1) = ts';
                    rpar(length(ts)+length(ws)+2,1) = reg;
                    rpar(length(ts)+length(ws)+3,1) = tal;
%                     rpar(length(ts)+length(ws)+4 : length(ts)+length(ws)+4 + lengrh(u)) = u;
                    
                    [zp0,ires] = state0(0,z0,zeros(length(z0),1),rpar'); % somente para EDO's
                    [tspan,zs,zps,outp]=dasslc('state1',[ti(ks):(ti(ks+1)-ti(ks))/20: ti(ks+1)],z0,zp0',rpar',rtol,atol);
                    z0 = zs(end,:)';
                    %--------------------------------------------------- 
%                 ts(end)
%                  pAD = myAD(u);
                    % ode...
%                  [tspan,zs]  = ode45(@(t,x)modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,...
%                      optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/2: ti(ks+1)], z0, optODE );
%                    [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,...
%uws,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks):(ti(ks+1)-ti(ks))/50:ti(ks+1)], z0, optODE );

                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                 t_out= [ t_out; tspan];
                 z0 = zs(end,:);
                 zss= [zss;zs];
        end
     else   % if est=0, full integration
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(ks,i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
        end
%            [tspan,zs]  = ode23( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
%             [tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [0 ts(ns+1)], z0, optODE );
               
                    rpar(1:length(ws),1) = ws(:,1);  % it is a vector....
                    rpar(length(ws)+1,1) = ks;
                    rpar(length(ws)+2:length(ws)+length(ts)+1,1) = ts';
                    rpar(length(ts)+length(ws)+2,1) = reg;
                    rpar(length(ts)+length(ws)+3,1) = tal;
                 
                    [zp0,ires] = state1(0,z0,zeros(length(z0),1),rpar'); % somente para EDO's
                    [tspan,zs,zps,outp]=dasslc('state1',[0:0.01:ti(end)],z0,zp0',rpar',rtol,atol);
                    z0 = zs(end,:)';
            
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end

            ko = length(tspan)+ko;
            t_out= [ t_out; tspan];
            z0 = zs(end,:);
            zss= [zss;zs];
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
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
 
%%
else
%% Forward state integration with Analytic (user specified) Jacobian...
  z0_emso = x0; % only for emso simulator
  zss= [];
  zdf = ones(n_s*ns,1);
  zdf = 0*zdf;
  z0 = [ x0' ; zdf ];

  if frozen==0
%     check = 2:length(ts);
    for i=1:n_c
        for j=2:ne(i)
            check(i,j) = j;
        end
    end
 end 
 
 
   %#####################################################################
 % DassL MODEL
 %#####################################################################
 dae_index = 0;
 Rz0 = [ x0 ]; 
 if reg==0
     ti = unique(abs(ts));  %organizing time steps...
     ko=0;
     
     
              % parameters rpar for dassl solver ...
               rpar(1,1) = n_c;
               rpar = [n_c length(ts) 0 0 ts  reg n_c n_s ts_frozen check  ns_a length(check)];
               rpar=rpar';
               
            for ks = 1:length(ti)-1

                
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i));   % calculating control variable u...
                end
            
                  % parameters rpar for dassl solver inside the loop...
                   rpar(3:length(u(i))+2,1) = u(:,1);  
                   rpar(length(u(i))+3,1) = ks;
                   
                  
                   [zp0,ires] = model_jac(0,z0,zeros(length(z0),1),rpar'); % only for ODEs...
                    yp0=zp0(1:n_s);
                    y0=z0(1:n_s);
                   %[tspan,zs]  = ode45( @(t,x)sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check), [ti(ks) ti(ks+1)], z0, optODE );
                   [t,zs1,zps,outp]=dasslc('model',[ts(ks) ts(ks+1)],y0,yp0',rpar',rtol,atol,index(dae_index+1,:),'Problem_A.dat','jac');  % only EADs with Jacobian
                   [t,zs2,zps,outp]=dasslc('model_jac',[ts(ks)  ts(ks+1)],z0,zp0',rpar',rtol,atol);
                   zs = [zs1 zs2(:,n_s+1:length(zp0))];
                  z0 = zs(end,:)';
                 
                 
                [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,uws,tal,...
                    est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks):(ti(ks+1)-ti(ks))/2: ti(ks+1)], z0, optODE );
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                t_out= [ t_out; tspan];
                z0 = zs(end,:)';                
                zss= [zss;zs];
            end

 
 end
  %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 for is = 1:ns
 df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 f = zss(:,:);
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
                 [tspan,zs]  = ode45(@(t,x)modelo_eq(t,x,u' ,ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,...
                     optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks):(ti(ks+1)-ti(ks))/2: ti(ks+1)], z0, optODE );
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
            [tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [0 ts(ns+1)], z0, optODE );
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
            t_out= [ t_out; tspan];
            z0 = zs(end,:);
            zss= [zss;zs];
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