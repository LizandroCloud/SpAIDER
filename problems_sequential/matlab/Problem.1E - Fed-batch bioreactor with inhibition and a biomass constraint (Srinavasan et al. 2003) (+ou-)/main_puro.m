 % Lizandro de Sousa Santos
 % Doutorado em Engenharia Química - Programa de Engenharia Química - Coppe - UFRJ.
 % Primeiro código embrionário com uso de wavelets para adaptação das
 % malhas...
 
 %Problema-exemplo Reator Isortemico de Schlegel....(1993) 
 
 %clear all;
 %clf;
 format long;
 close all;
 % Opcoes do NLP
  optODE = odeset( 'RelTol', 1e-10, 'AbsTol', 1e-10 );
  optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 100000);

 % Horizonte de  tempo....
 t0 = 0; % Tempo Inicial
 tf = 250; % Tempo Final
 x0 = [ 0.72; 0.05; 1]; %Seis estimativas iniciais
 fob = -0.4; %Valor inicial da funcao-objetivo (para plotagem);
 ti = 0; %(instante inciial para plotagem);
for is = 7:7 % Indexador para variar o numero de estagios...
   
  %if is==1
    ns = 2^(is);
    wL = 0*ones(ns,1);  % limite minimo da u
    wU = 0.001*ones(ns,1);  % limite maximo da u
    w0 = ones(ns,1)*0.0005;  % estimativa inicial da u...
    ts = [t0:(tf-t0)/ns:tf]; % Equiparticão padrão...
 %end
 
 iu=is+2;  % é um contador para  
 % Chamando o FminCon...
 tic
 [wopt, erropt, iout] = fmincon( @(ws)obj(x0,ns,ts,ws,optODE), w0, [], [], [],...
 [], wL, wU, @(ws)ctr(x0,ns,ts,ws,optODE), optNLP);
 dispopt( x0, ns, ts, wopt, optODE,iu );
 toc
 t = [];

 
 fob(is+1) =  erropt;   % evolução da função objetivo...
 ti(is+1) = is;
 %if ( norm( fob(is+1)- fob(is) )/fob(is))<0.05
 %    break;
 %end

end
 


 
  
  

 
 
 