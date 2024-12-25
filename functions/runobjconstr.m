function [wopt, erropt, iout,output,lambda,grad,hessian] =...
    runobjconstr(w00, wLL, wUU, nonl_constraint , optNLP,x0,nss,tss,ns,Nxo,nUt,w_frozen,check,...
    ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,in_SPAIDER)

global Jx;
global dJx;
global ceqx;
global cx;
global dceqx;
global dcx;
global flast;
if nargin == 1 % No options supplied
    optNLP = [];
end

wsLast = []; % Last place computeall was called
% myf = []; % Use for objective at xLast
% myc = []; % Use for nonlinear inequality constraint
% myceq = []; % Use for nonlinear equality constraint
dflast=[];
flast=[];
toutL=[];
tuL=[];
yL=[];
d2last = []; 
nL = [];


switch (solver)  % for each solver
      case 'mbo'
     solve = '@ode15i';
     objective = @(ws)objfunc_mbo(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solve,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_mbo(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solve,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'emso'
     solver.func = @emso;
     objective = @(ws)objfunc_emso(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_emso(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode45'  
         
      solver = '@ode45';
      objective = @(ws)objfunc_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
      case 'ode15i'  
         
      solver.func = @ode15i;
      objective = @(ws)objfunc_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode23s'   
     solver.func = @ode23s;    
     objective = @(ws)objfunc_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode15s'
     solver= @ode15s;      
     objective = @(ws)objfunc_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'ode23t'
     solver.func = @ode23t;       
     objective = @(ws)objfunc_odex(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_ode(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     case 'dasslc' 
     solver.func = @dasslc;
     objective = @(ws)objfunc_dasslc(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER);
     nonl_constraint = @(ws)constr_dasslc(x0,nss,tss,ws,ns,Nxo,nUt,w_frozen,check,...
         ns_a,frozen,is,is0,w0,ts_frozen, solver.func,optODE,typ,reg,est,wLL, wUU,in_SPAIDER); 
end

[wopt, erropt, iout,output,lambda,grad,hessian] = ....
      fmincon(  objective, w00, [], [], [], [], wLL, wUU, nonl_constraint , optNLP);   
% nvars = (ns(1)+ns(2))*nUt;
% options = optimoptions('particleswarm','SwarmSize',20,'HybridFcn',@fmincon, 'Display', 'iter', 'UseParallel', false, 'FunctionTolerance', 1e-3, 'MaxStallIterations', 10);
% grad=[]; hessian=[]; lambda=[];
% [wopt, erropt, iout, output, lambda] = particleswarm(objective,nvars,wLL,wUU,options);
% options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
% [wopt, erropt, iout, output] =  patternsearch(objective,w00,[],[],[],[],wLL, wUU, nonl_constraint,options);

 function [ J, dJ ] = objfunc_odex( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER)

 wsN=ws;
 k=0; ki=1;  t_ult=[]; t_out=[]; yout=[]; n_elem=[];
 if ~isequal(ws,wsLast)  % Check if computation is necessary
     if nargout == 1
         [f, t_ult, t_out, yout,n_elem] = function_ode( x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
         wsLast = ws;
         flast = f;
         tuL=t_ult;
         toutL = t_out;
         yL = yout;
         d2last = 1; % checking if df was atualized...
         nL = n_elem;
     else
         [f,t_ult,t_out,yout,n_elem,df] = function_ode(  x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
        wsLast = ws;
        flast = f;
        dflast = df;
        tuL=t_ult;
        toutL = t_out;
        yL = yout;
        d2last = 2; % checking if df was atualized...
        nL = n_elem;
     end
 end
     if nargout == 1
        f  = flast;
        ws = wsLast;
        t_ult=tuL;
        t_out=toutL;
        yout=yL;
        [ c, ceq] = constr_odex( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER  );
        [J]= obf1(f,t_ult,t_out,ws,yout) ;
        for i=1:length(c)
            sd(i) = max(c(i),0);
        end
        for j=1:length(ceq)
            sh(i) = max(ceq(j),0);
        end
        % J=J+norm(sd)*10 + norm(sh)*10;
     else
        f  = flast;
        ws = wsLast;
        t_ult=tuL;
        t_out=toutL;
        yout=yL;
            if d2last==2
                df = dflast;
            else
                [f,t_ult,t_out,yout,n_elem,df] = function_ode(  x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER); 
            end
        [J dJ]= obf2(f,t_ult,t_out,df,wsN,yout);
         dflast = df;
         d2last = 2;
         nL = n_elem;
     end

 end


 function [ c, ceq, dc, dceq ] = constr_odex( x0, ns, ts, ws, ne, n_s,n_c,w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,wLL, wUU,optSPAIDER  )
 c = []; ceq = []; dc = []; dceq = [];  t_ult=[]; t_out=[]; yout=[]; n_elem=[];

 wsN=ws;
 k=0; ki=1;
 
 if typ==1 % if terminal free
     wsN(n_c+1)=ws(n_c+1)*(1000-100)+100;
 end

  if ~isequal(ws,wsLast)  % Check if computation is necessary
     if nargout == 2 
         [f, t_ult, t_out, yout,n_elem] = function_ode( x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
         wsLast = ws;
         flast = f;
         tuL=t_ult;
         toutL = t_out;
         yL = yout;
         d2last = 1; % checking if df was atualized...
         nL = n_elem;
     else
         [f,t_ult,t_out,yout,n_elem,df] = function_ode(  x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER);
        wsLast = ws;
        flast = f;
        dflast = df;
        tuL=t_ult;
        toutL = t_out;
        yL = yout;
        d2last = 2; % checking if df was atualized...
        nL = n_elem;
     end
 end
     if nargout == 2
        f  = flast;
        ws = wsLast;
        t_ult=tuL;
        t_out=toutL;
        yout=yL;
        n_elem = nL;
        [ceq c] = cr1(f, t_ult,t_out,ws,yout,n_elem,ne);  % Modificação UOTR
     else
        f  = flast;
        ws = wsLast;
        t_ult=tuL;
        t_out=toutL;
        yout=yL;
        n_elem = nL;
            if d2last==2
                df = dflast;
            else
                [f,t_ult,t_out,yout,n_elem,df] = function_ode(  x0, ns, ts, ws, ne,n_s,n_c,...
             w_froz,check,ns_a,frozen,is,is0,w0,ts_frozen,solver,optODE,typ,reg,est,optSPAIDER); 
            end
        [ceq, c, dceq, dc] = cr2(f, t_ult,t_out,df,ws,yout,n_elem, ne);  % Modificação UOTR
         dflast = df;
         d2last = 2;
     end
 
 end


end