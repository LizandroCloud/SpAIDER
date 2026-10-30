 clear all;
 clf;
 format long;

 % Options for ODE & NLP Solvers
 optODE = odeset( 'RelTol', 1e-9, 'AbsTol', 1e-9 );
 optNLP = optimset( 'LargeScale', 'on', 'GradObj', 'on', 'GradConstr', 'on',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-9),...
 'TolFun', 10^(-9), 'TolCon', 10^(-9), 'MaxFunEval', 100);

 % Time Horizon and Initial State
 t0 = 0; % Initial time
 tf = 2; % Final time
 x0 = [ 1; 1; 0 ]; % Initial state

 ns = 2;
 for is = 1:3
     
 ns = 2^ns; % Numero de estagios: ns = 2, 4, 8, 16, and 32
 
 if is==1
    wL = -10*ones(ns,1);
    wU = 10*ones(ns,1);
    w0 = zeros(ns,1);
    ts = [t0:(tf-t0)/ns:tf]; % Time stages (equipartition)
 else
    w0 = waverec(C,L,'db1');
    ts(1) = 0.0;
    j =1;
    for i=1:(ns-1)
        if w0(i+1)==w0(i)
        else
        ts(j+1)=t0+((tf-t0)/ns)*i;
        ww0(j)=w0(i);
        j = j+1;
        end
    end
 ts(j+1)=tf;
 ww0(j)=w0(int32(ns));
 ns = j; % atualizo ns novo irregular.
 wL = -10*ones(ns,1);
 wU = 10*ones(ns,1);
 w0=ww0; % atualizo novo perfil irregular.
 end


 % CVP-based Optimization
 [wopt, erropt, iout] = fmincon( @(ws)obj(x0,ns,ts,ws,optODE), w0, [], [], [],...
 [], wL, wU, @(ws)ctr(x0,ns,ts,ws,optODE), optNLP);
 dispopt( x0, ns, ts, wopt, optODE );
 
 for i=1:(length(wopt)-1)
     Cjj(2*i)=(wopt(i)+(wopt(i+1)))/2.0;
 end 
 for i=1:2:(length(wopt)*2)
     Cjj(i)=wopt((i+1)/2);
 end 
     Cjj(2*ns) = wopt(ns)+ (wopt(ns-1)+(wopt(ns)))/2.0;
 [C,L] = wavedec(Cjj,int32(log2(ns)),'db1');
  
 for i=1:(length(C))
 if sqrt(C(i)^2.0)<=1.0;
     C(i)=0.0;
 end
 end
  
 ns =  log2(ns) + 1;
 
 end
 
 
 