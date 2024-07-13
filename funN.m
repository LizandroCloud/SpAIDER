function [ f, t_out, df ] = fun( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ )
 
global reg;
global est;
global DEBUG;  % for dassl solver...
global t_ult;

%% Parameters
DEBUG  = 1; 
dae_index=0;
t0 = 0;
ks=1;
reg = 0; % regularization...
fix =0;  % fix solver...
est = 1;  % integration by steps...
ts = roundn (ts, -5);
koz=1;  %reorganizing time steps...
k=1;
ko=0;
ne_a=ne;
rtol= 1e-10;       % relative tolerance for dassl
atol= 1e-10;       % absolute tolerance for dassl


if typ==1
      tf=ws(end);
    if isempty(t_ult)
        t_ult=tf;
    end
  
    for i=1:length(ts)
        ts(i)=ts(i)*ws(end)/t_ult;
    end
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
if nargout == 2 % Forward state integration with numerical jacobian...

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
 if strcmp(solver,'fix_runge')   % If its a fixed solver....
     if est==1  % with integration by steps...
         for ks = 1:ns
             [zs] = ode5( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks):0.5:ts(ks+1)], [z0]);
             z0 = zs(end,:)';
         end
     else  % complete integration...
             [zs] = ode5( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks):0.3:ts(ns+1)], [z0]);
     end
 elseif    strcmp(solver,'ode')   
     if est==1  % with integration by steps...
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
%                 ts(end)
                 [tspan,zs]  = ode45( @(t,x) modelo_eq(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check), [ti(ks) ti(ks+1)], z0, optODE );
                 t_out= [ t_out; tspan];
                 z0 = zs(end,:);
                 zss= [zss;zs];
        end
     else   
           [tspan,zs]  = ode45( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
     end
 elseif  strcmp(solver,'dassl') 
     index = zeros(1,n_s);
     if est==1  % with integration by steps...
                
               ti = unique(abs(ts));  %organizing time steps...
          
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
                   
                  
                   [zp0,ires] = model(0,z0,zeros(length(z0),1),rpar'); % only for ODEs...
                    yp0=zp0(1:n_s);
                    y0=z0(1:n_s);
                   %[tspan,zs]  = ode45( @(t,x)sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check), [ti(ks) ti(ks+1)], z0, optODE );
%                    [t,zs,zps,outp]=dasslc('model',[ts(ks):0.001:ts(ks+1)],y0,yp0',rpar',rtol,atol,index(dae_index+1,:),'Problem_A.dat','jac');  % only EADs with Jacobian
                   [t,zs,zps,outp]=dasslc('model',[ts(ks) ts(ks+1)],y0,yp0',rpar',rtol,atol); 
%                    zs = [zs1 zs2(:,n_s+1:length(zp0))];
                   z0 = zs(end,:)';
                   zss= [zss;zs];
            end
         
         
     else   
           [tspan,zs]  = ode45( @(t,x)model(t,x,ws,ks,ts,reg), [ts(ks) ts(ns+1)], z0, optODE  );
     end
 elseif  strcmp(solver,'hysys') 
      if est==1  % with integration by steps...
        ti = unique(abs(ts));  %organizing time steps...
        for ks = 1:length(ti)-1
                for i = 1:length(wx(:,1)) 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
%                 ts(end)
                 [tspan,zs]  = hysys( [ti(ks) ti(ks+1)], u, z0 );
                 t_out= [ t_out; tspan];
                 z0 = zs(end,:);
                 zss= [zss;zs];
        end
     else   
		    ti = unique(abs(ts));  %organizing time steps...
           [tspan,zs]  = hysys( [ti(1) ti(end)], u, z0 );
     end
 
 end
      
 % Final values
 t_out;
 f = zss(:,:);
 t_ult=t_out(end);
%%
else
%% Forward state integration with Analytic (user specified) Jacobian...
  zss= [];
  zdf = ones(n_s*ns,1);
  zdf = 0*zdf;
  z0 = [ x0'; zdf ];

  if frozen==0
%     check = 2:length(ts);
    for i=1:n_c
        for j=2:ne(i)
            check(i,j) = j;
        end
    end
 end 
 
 if  strcmp(solver,'ode')
 
 Rz0 = [ x0 ];
 if reg==0
     ti = unique(abs(ts));  %organizing time steps...
            for ks = 1:length(ti)-1
                
             
                for i = 1:n_c 
                  u(i) = u_reg1(0.5*(ti(ks)+ ti(ks+1)),te(i,:),wx(i,:),ks,ne(i)); 
                end
                
                for j=1:n_c
                     for i = 1:ne(j) 
                        tx = [0.5*(ti(ks)+ ti(ks+1)) ts_frozen(j,2:ne(j))];
                        tx = sort(tx);
                        posi = find(tx == 0.5*(ti(ks)+ ti(ks+1)),1,'last');
                
                
                        if ((posi)== i)
                            uws(j,i) = 1;
                            else   % do contrario...
                            uws(j,i) = 0;
                        end
                     end
                end
                
                [tspan,zs]  = ode45( @(t,x) modelo_eq_sens(t,x,u',ks,ts,reg,ns,ne,te,n_s,n_c,ne_a,ts_frozen,check,uws), [ti(ks) ti(ks+1)], z0, optODE );
                 t_out= [ t_out; tspan];
                z0 = zs(end,:)';                
                zss= [zss;zs];
            end
 else
    %    for ks = 1:ns
    [tspan,zs] = ode45( @(t,x)sens(t,x,ws,ks,ts,reg,ne), [0 0.2], z0, optODE );
    %    z0 = zs(end,:)';
    %  end
 end
 
 else  % if solver is dasslc
 index = zeros(1,n_s);
 if reg==0
     ti = unique(abs(ts));  %organizing time steps...
     
     
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
            end
 else
    %    for ks = 1:ns
    [tspan,zs] = ode45( @(t,x)sens(t,x,ws,ks,ts,reg,ne), [0 0.2], z0, optODE );
    %    z0 = zs(end,:)';
    %  end
 end
 
 end
 
  t_ult=t_out(end);
 % keeping the objective function value
 f = zss(:,:);


 
 for is = 1:ns
 df(:,1:n_s,is) = zss(:,is*n_s+1:is*n_s+n_s);
 end
 end
end