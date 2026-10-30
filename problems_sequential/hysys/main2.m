% Problema HSIAO


 % Opcoes do NLP
  optODE = odeset( 'RelTol', 1e-7, 'AbsTol', 1e-7 );
  optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-6),...
 'TolFun', 10^(-6), 'TolCon', 10^(-6), 'MaxFunEval', 100000);
  close all;
%  'rigrsure' use principle of Stein's Unbiased Risk.  ...Mais conservativa
%  'heursure' is an heuristic variant of the first option.
%  'sqtwolog' for universal threshold sqrt(2*log(.)).
%  'minimaxi' for minimax thresholding.  ...Mais conservativa
  
 opC = 1;  % 1 - Wav ; 2 -> ADE
 % Horizonte de tempo....
 t0 = 0; % Initial time
 tf = 2; % Final time
 x0 = [ 1; 1; 0 ]; % Initial state

 wname = 'haar';  % tipo de base wavelets.
 tptr = 'sqtwolog';
 sorh = 's';
 scal = 'mln';
 is0 = 3;
 isF = 3;
 thrs = [];
 w=[];
 w(1:2^isF,1:isF)=0;
 tx(1:2^isF,1:isF)=0;
for is = is0:isF

% if is<isF
%  optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
%  'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-5),...
%  'TolFun', 10^(-5), 'TolCon', 10^(-5), 'MaxFunEval', 100000);
% else
%   optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
%  'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
%  'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 100000);
% end
% Numero de estagios: ns = 2, 4, 8, 16, and 32 ... 2^is
Cgrad=0;
grad=0;
Dgrad=0;
Pgrad=0;
if is==is0
   ns = (2^is);
     wL = -10*ones(ns,1);
    wU = 10*ones(ns,1);
    w0 = ones(ns,1);
    ts = [t0:(tf-t0)/ns:tf]; % Time stages (equipartition)
    point = [1:ns]*0;
end
 iu=is+2;
 % Chamando o FminCon...
 tic % contagem do tempo...
 [wopt, erropt, iout] = fmincon( @(ws)obj(x0,ns,ts,ws,optODE), w0, [], [], [],...
 [], wL, wU, @(ws)ctr(x0,ns,ts,ws,optODE), optNLP);
 toc
 
 % armazenamento da variavel de controle em varios estagios do Awav...
 
 w(1:length(wopt),is-2)=wopt;
 tx(1:length(ts),is-2)=ts;
 
 normGrad=norm(wopt(2:(length(wopt)-1)))/(length(wopt)-2);
 for ih=2:(length(wopt))
    Dgrad(ih-1)= abs(wopt(ih))-abs(wopt(ih-1));
 end
 il=1;
 
 if opC==1
 valK = 0.2;
 else
  valK = 0.0;
 end
 for iw=1:length(Dgrad)
%      Dgsin=Dgrad(iw)*Dgrad(iw-1);
     if abs(Dgrad(iw))>(valK)*normGrad
         Cgrad(il)= iw+1;
         il=il+1;
     end
 end
 
 % pontos prospectivos...
 Pgrad(1)= Cgrad(1)+1;
%   Pgrad(2)= Cgrad(1)+2;
 Pgrad(2)= Cgrad(length(Cgrad))-2;
%   Pgrad(4)= Cgrad(2);
 Pgrad=Cgrad;
 
 dispopt( x0, ns, ts, wopt, optODE,iu ); % plotando...
 
 if opC==1 % Analise wavelets
 maxlev = WMAXLEV(length(wopt),wname);
 t = [];
 figure(iu+530);
% for i=2:ns+1
% u(i) = wopt(i-1);
% end
% plot(ts,u,'o');
 if is>is0
 wopt = wopt';
 % Colocando no dominio wavelets e removendo pontos...
 [C,L] = wavedec(wopt,maxlev,wname);  % C sao os detalhes...
 tpq = wpdec(wopt,maxlev,wname);
 for i=1:length(C)
  ci(1,i)=C(i);
 end
 for i=1:length(L)
  li(1,i)=L(i);
 end

   
 figure(iu+1000+530); % plotando os coeficientes...
 cfd = zeros(maxlev,length(wopt)); 
 for k = 1:maxlev
    d = detcoef(ci,li,k); 
    d = d(ones(1,2^k),:); 
    cfd(k,:) = wkeep( d(:)',length(wopt) ); 
 end 

 cfd = cfd(:); 
 I = find(abs(cfd)<sqrt(eps)); 
 cfd(I)=zeros(size(I)); 
 cfd = reshape(cfd,maxlev,length(wopt));

% Plotando os coeficientes. 
subplot(312), colormap(1-gray(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
title('Coeficientes da transformada wavelets') 
ylabel('Nivel')

j=1;
%  Nc = norm(C);  % Norma dos coeficientes...
%     for i=L(1)+1:length(C)
%       if abs(C(i))<1e-4;  %critério para eliminação de pontos...
%         C(i)=0.0;
%       end
%     end
[w0,C,L] = wden(wopt,tptr,sorh,scal,maxlev,wname);
 thr = thselect(wopt,tptr);  % valor do threshold...
thrs(is)=thr;

%[C,L] = wavedec(w0,maxlev,wname);
 for i=1:length(C)
  ci(1,i)=C(i);
 end
 for i=1:length(L)
  li(1,i)=L(i);
 end
%tpq = wpdec(w0,maxlev,wname);
%  figure(iu+4000);
%plot(tpq);
figure(iu+2000+530);
w0 = waverec(C,L,wname); % transformada wavelets inversa...
cfd = zeros(maxlev,length(wopt)); 
 
for k = 1:maxlev  % Organizando o gráfico com os coeficientes de cada nivel
  d = detcoef(ci,li,k); 
  d = d(ones(1,2^k),:); 
  cfd(k,:) = wkeep(d(:)',length(wopt)); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(wopt));

% Plot discrete coefficients. 
subplot(312), colormap(1-gray(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
title('Coeficientes da transformada wavelets sem ruido') 
ylabel('Nivel')
 % dividindo estagios...
    
 else  % se for menor que is0
     w0=wopt;
     w0ns=wopt;
end
 else % sem wavelets
     if is>is0
        for i = 1:length(wopt)
            wopt(i) = wopt(i) + rand(1)*1e-4;  % garantir pequena diferença entre valores
            w0=wopt;
        end;
     else
         w0=wopt;
         w0ns=wopt;
     end
 end
 jj=1;

for j=1:2:ns-1  % dããã...  Da pra melhorar essa lógica...(vou fazê-lo)
        for hs=1:length(Pgrad)  % verificando se o Grad atual é o selecionado...
                delta(hs)=((w0(j))-w0(Pgrad(hs)));
                if delta(hs)==0
                    ind=1;
                end
        end
       if w0(j)==w0(j+1)  % depois da transformação...   
            w0ns(jj)=w0(j);      %nao vou dividir nesse caso...
            w0ns(jj+1)=w0(j+1); 
            ts(jj+1)=t0+((tf-t0)/(2*ns))*(2*j);  %
            ts(jj+2)=t0+((tf-t0)/(2*ns))*( (2*j+2) );  %
            point(jj)=1; point(jj+1)=1;
            jj=jj+2;
       else
            ind=0;
            if length(w0)<2
            for hs=1:length(Pgrad)  % verificando se o Grad atual é o selecionado...
                delta(hs)=((w0(j))-w0(Pgrad(hs)));
                if delta(hs)==0
                    ind=1;
                end
            end
            else
                for hs=1:length(Pgrad)  % verificando se o Grad atual é o selecionado...
                    delta(hs)=((w0(j+1))-w0(Pgrad(hs)));
                    if delta(hs)==0
                        ind=1;
                    end
                end
                for hs=1:length(Pgrad)  % verificando se o Grad atual é o selecionado...
                delta(hs)=((w0(j))-w0(Pgrad(hs)));
                if delta(hs)==0
                    ind=1;
                end
                end
             end
                if (ind==1)
                     w0ns(jj)=w0(j);   w0ns(jj+1)=w0(j);    %aqui divido...
                     w0ns(jj+2)=w0(j+1); w0ns(jj+3)=w0(j+1);
                     ts(jj+1)=t0+((tf-t0)/(2*ns))*(2*(j-1)+1); ts(jj+2)=t0+((tf-t0)/(2*ns))*(2*(j-1)+2);
                     ts(jj+3)=t0+((tf-t0)/(2*ns))*(2*(j-1)+3); ts(jj+4)=t0+((tf-t0)/(2*ns))*(2*(j-1)+4);
                     jj=jj+4;
                else
                     w0ns(jj)=w0(j);      %nao vou dividir nesse caso...
                     w0ns(jj+1)=w0(j+1); 
                     ts(jj+1)=t0+((tf-t0)/(2*ns))*(2*j);  %
                     ts(jj+2)=t0+((tf-t0)/(2*ns))*( (2*j+2) );  %
                     point(jj)=1; point(jj+1)=1;
                     jj=jj+2;
                end
       end
end      
 ns = length(w0ns);
 wL = -10*ones(ns,1);
 wU = 10*ones(ns,1);
 w0 = zeros(ns,1);
 %wL = -10*ones(ns,1);
 %wU = 10*ones(ns,1);
 
 li=0;
 ci=0;
 end




 
  
  

 
 
 