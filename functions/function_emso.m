function [ f, t_ult, t_out, y, ti, df ] = function_emso( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER )
 
% global reg;
% global est;
global DEBUG;  % for dassl solver...
% global t_ult;

%% Parameters
DEBUG  = 1; 
dae_index=0;
t0 = 0;
ks=1;
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
ne_ms= length( unique(abs(ts)) ) -2;
% ne_ms=(dim(2)-n_c*sum(ne))/(n_s) -0; % numer of elementos of multiple shooting
if strcmp(optSPAIDER.problem.method,'multiple shooting')  % keeping the initial conditions for multiple shooting
    for j=1:n_s
    x0_ms(1:ne_ms,j) = ws(sum(ne) + (j-1)*ne_ms + 1 : sum(ne) + (j)*ne_ms);
    end
end


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
 
     %#####################################################################
     % EMSO MODEL
     %#####################################################################
     
     ko=0;
     if est==0  % if the integration is from to to tf...
        ti = unique(abs(ts));  %organizing time steps (using unique to classify)...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1))  % length of each control vector...
                  u(ks,i) = u_reg1(1*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); % u_reg1 is a function heaviside regularized...
                end
        end
%                 ts(end)
%                  [tspan,zs]  = ode23s( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est), [ti(ks) ti(ks+1)], z0, optODE );
%                  tic
               %coment    [tspan,hh,zs] = sim('emso',[0  ti(end)],[],[ti(2:end)' u ones(length(ti)-1,1)*z0]);
                   u=[u ;u(end)];
                   pause(5)
                   try
                   [tspan,hh,zs] = sim('fcn2',[0:10:  ti(end)],[],[ti(1:end)' [u,u]]);
                   catch
                       disp('Erro na simulação')
                       set_param('fcn2B','TryForcingSFcnDF','on')
                       [tspan,hh,zs] = sim('fcn2B',[0:10:  ti(end)],[],[ti(1:end)' [u,u]]); 
                       %set_param('fcn2','TryForcingSFcnDF','on')
                       %zs = zs';
                   end
                 % [tspan,hh,zs] = sim('fcn',[0  ti(end)],[],[ti(2:end)' u]);
 %                   [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0  ]);
%                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0 ones(length(ti)-1,1)*tal ]);
%                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u  ]);
%                  toc
%                 [373.730949988220,0.264044804775476,0.206912759316386,0.263994530507903,0.00108039002044736,0.263967515379788,0.000352522413560523,0.181475950519473,7.31818291362804e-06,0.817492055838849,0.000672153045204316;373.727409054339,0.264044596941463,0.206918634662591,0.263988393532297,0.00108044829998764,0.263967926563661,0.000352532810600101,0.181474052953428,7.31923386201321e-06,0.817493928942709,0.000672166059400157;]
                
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                 t_out= [ t_out; tspan];
%                  z0 = zs(end,:);
                 zss= [zss;zs];                     
       % end
     else   % if est=1, integration by steps (or elements)
         
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
                for (ks = 1:tix)
                        u={}; nex={}; i=0; ugx={}; tex={};
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
                         [tspan,hh,zs]= parsim('emso',[0 :0.5: 1],[  tf-to]', [u{ks,:} dtt pzo]); % original
%                           [tspan,hh,zs] = sim('emso',[ti(1) :0.5: 1],[],[(ti(ks+1)-ti(ks))' u (ti(ks+1)-ti(ks)) ones(1,1)*z0 ]);
%                        [tspan,hh,zs] = sim('emso',[ti(1) :0.5: 1],[],[(ti(ks+1)-ti(ks))' u (ti(ks+1)-ti(ks)) ones(1,1)*z0 ]);
%                [tspan,hh,zs] = sim('emso',[ti(ks) :( (ti(ks+1)-ti(ks))/10 ): ti(ks+1)],[],[ti(ks+1)' u ones(1,1)*z0]);
                 %[tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks) ti(ks+1)], z0, optODE );
                         tspan = tspan*dtt;
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
%      toc
     zss=[]; t_out=[];
     for kx=1:tix
         zsa =  (zsc{1,kx});
         zss = [zss;zsa];
         tspan = (tsc{kx});
         t_out= [ t_out; (kx-1)*tspan(end)+tspan];
     end

     for ivf = 1:tix
                for i = 1:length(wx(:,1))  % length of each control vector...
                  ui(ivf,i) = u_reg1(1*(ti(ivf)+ ti(ivf+1)),te(i,:),wx(i,:),1,ne(i)); % u_reg1 is a function heaviside regularized...
                end
         [yout] = alg_eq(zss(ivf,:),ui);
         y(ivf,:) = yout;
     end
      
        else % single shooting
              ti = unique(abs(ts));  %organizing time steps...   
            for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i));  % calculating the value of u in time domain (ti)
                end
                   dtt = (ti(ks+1)-ti(ks));
%                  [tspan,hh,zs] = sim('emso',[ti(1) :( (ti(ks+1)-ti(ks))/10): ti(ks+1)],[],[ti(ks+1)' u ]);
%                    [tspan,hh,zs] = sim('emso',[ti(1) :1/10: (ti(ks+1)-ti(ks))],[],[(ti(ks+1)-ti(ks))' u  ones(1,1)*z0 ]);
                    [tspan,hh,zs] = sim('emso',[ti(1) :0.5: 1],[],[(ti(ks+1)-ti(ks)) u dtt ones(1,1)*z0]);
                 %[tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks) ti(ks+1)], z0, optODE );
                  tspan = tspan*(dtt);
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                 t_out= [ t_out; (ks-1)*tspan(end)+tspan];
                 z0 = zs(end,:);
                 zss= [zss;zs(1:end-1,:)];
            end


        end
        
     end    
     
     %#####################################################################
 
%  toc
 %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 f = zss(:,:);
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
 
%     ti = unique(abs(ts));  %organizing time steps (using unique to classify)...
%         for ks = 1:length(ti)-1
%                 for i = 1:length(wx(:,1))  % length of each control vector...
%                   u(ks,i) = u_reg1(1*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); % u_reg1 is a function heaviside regularized...
%                 end
%         end
% %                 ts(end)
% %                  [tspan,zs]  = ode23s( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est), [ti(ks) ti(ks+1)], z0, optODE );
% %                  tic
%                    a=1;
%                    p = rand(5,1);
%                    p = [0    1 5.0000    1.0000   41.6929]';
%                    p = [0    1 ]';
%                    pAD = myAD(p);
% %                    pAD = myAD([u ;(ones(length(ti)-1,1)*z0)']);
%                    
%                    [tspan,hh,zs] = sim('emso',[0  ti(1)],[],[ti(1:2)' pAD]);
%                    
%                 
%  %                   [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0  ]);
% %                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0 ones(length(ti)-1,1)*tal ]);
% %                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u  ]);
% %                  toc
% %                 [373.730949988220,0.264044804775476,0.206912759316386,0.263994530507903,0.00108039002044736,0.263967515379788,0.000352522413560523,0.181475950519473,7.31818291362804e-06,0.817492055838849,0.000672153045204316;373.727409054339,0.264044596941463,0.206918634662591,0.263988393532297,0.00108044829998764,0.263967926563661,0.000352532810600101,0.181474052953428,7.31923386201321e-06,0.817493928942709,0.000672166059400157;]
%                 
%                  for ivf = 1+ko:length(tspan)+ko
%                      [yout] = alg_eq(zs(ivf-ko,:),u);
%                      y(ivf,:) = yout;
%                  end
%                  ko = length(tspan)+ko;
%                  t_out= [ t_out; tspan];
% %                  z0 = zs(end,:);
%                  zss= [zss;zs];                     
%        % end  
%  
%  
 
 
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
 % EMSO MODEL
 %#####################################################################
 ko=0;
     if est==0  % if the intergration is from to to tf...
        ti = unique(abs(ts));  %organizing time steps (using unique to classify)...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1))  % length of each control vector...
                  u(ks,i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); % u_reg1 is a function heaviside regularized...
                end
        end
%                 ts(end)
%                  [tspan,zs]  = ode23s( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est), [ti(ks) ti(ks+1)], z0, optODE );
%                  tic
                   [tspan,hh,zs] = sim('emso',[0 0.001 ti(end)],[],[ti(2:end)' u ones(length(ti)-1,1)*z0_emso]);
 %                   [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0  ]);
%                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u ones(length(ti)-1,1)*z0 ones(length(ti)-1,1)*tal ]);
%                     [tspan,hh,zs] = sim('emso',[0:0.01: ti(end)],[],[ti(1:end-1)' u  ]);
%                  toc
%                 [373.730949988220,0.264044804775476,0.206912759316386,0.263994530507903,0.00108039002044736,0.263967515379788,0.000352522413560523,0.181475950519473,7.31818291362804e-06,0.817492055838849,0.000672153045204316;373.727409054339,0.264044596941463,0.206918634662591,0.263988393532297,0.00108044829998764,0.263967926563661,0.000352532810600101,0.181474052953428,7.31923386201321e-06,0.817493928942709,0.000672166059400157;]
                
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                 t_out= [ t_out; tspan];
                 z0 = zs(end,:);
                 zss= [zss;zs];                     
       % end
     else   % if est=0, full integration
%          
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i));  % calculating the value of u in time domain (ti)
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
                
                % dae model
                [tspan,hh,zse] = sim('emso',[ti(ks) :( (ti(ks+1)-ti(ks))/10 ): ti(ks+1)],[],[ti(ks+1)' u 1 ones(1,1)*z0_emso]);
                % only sensibility equations
                [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,uws,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.tref), [ti(ks) :( (ti(ks+1)-ti(ks))/10 ):  ti(ks+1)], z0, optODE );
                zs(:,1:n_s) = zse(:,1:n_s);
                 %[tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,tal,est,optSPAIDER.prob.uref,optSPAIDER.prob.tref,optSPAIDER.prob.nref), [ti(ks) ti(ks+1)], z0, optODE );
                 
                 for ivf = 1+ko:length(tspan)+ko
                     [yout] = alg_eq(zs(ivf-ko,:),u);
                     y(ivf,:) = yout;
                 end
                 ko = length(tspan)+ko;
                 t_out= [ t_out; tspan];
                 z0_emso = zse(end,:); % for emso model
                 z0 = zs(end,:);
                 zss= [zss;zs];
        end
     end
     %#####################################################################
     

     
     
     
     
     
     
   %+++++++++++++++++++++++++++++++++++++++++++++++
 % Final values for calculating objective function
 for is = 1:ns
 df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 f = zss(:,:);
 t_ult=t_out(end)*tal;
 %+++++++++++++++++++++++++++++++++++++++++++++++
end