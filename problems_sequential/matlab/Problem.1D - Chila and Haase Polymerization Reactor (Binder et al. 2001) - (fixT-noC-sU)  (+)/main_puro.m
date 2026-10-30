 % Lizandro de Sousa Santos
 % Doutorado em Engenharia Química - Programa de Engenharia Química - Coppe - UFRJ.
 % Primeiro código embrionário com uso de wavelets para adaptação das
 % malhas...
 
 %Problema-exemplo Reator Polimerizacao Chylla e Haase....(1993) 
 
 % Opcoes do NLP
optODE = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8 );
  optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 10000);
  close all;
  
 % Horizonte de tempo....
 t0 = 0; % Tempo Inicial
 tf = 18000; % Tempo Final
 x0 = [ 0; 11.227; 325; 0 ]; %Condição inicial: 1-mp 2-mn 3-T 4-Fobj
 fob = 6.5;  % Valor inicial da funcao-objetivo (para plotagem);
 ti = 0; %(para plotagem);
 
for is =2:2
    
    ns = 2^(is);
    wL = 270*ones(ns,1);  % limite minimo da u
    wU = 400*ones(ns,1);  % limite maximo da u
    x = t0:(((tf-t0)/(ns-1))):tf;
    y = 345*ones(ns,1);
    xx = linspace(t0,tf,8);
    cs = spline(x,y);
    yy = ppval(cs,xx);
    ts = [t0:(tf-t0)/ns:tf]; 
    w0 = yy;
  iu=is+3;
 tic
 [wopt, erropt, iout,output,lambda,grad,hessian] = fmincon( @(ws)obj(x0,ns,ts,ws,optODE), w0, [], [], [],...
 [], wL, wU, @(ws)ctr(x0,ns,ts,ws,optODE), optNLP);
 toc
 dispopt( x0, ns, ts, wopt, optODE,iu );
 t = [];
 

 
fob(is+1) =  erropt;   % evolução da função objetivo...
  ti(is+1) = is;
 if ( norm( fob(is+1)- fob(is) )/fob(is))<0.05
     break;
 end

end


 
  
  

 
 
 