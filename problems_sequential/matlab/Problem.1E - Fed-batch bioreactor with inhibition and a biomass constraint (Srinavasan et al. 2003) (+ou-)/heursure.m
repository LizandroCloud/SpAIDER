 % Lizandro de Sousa Santos
 % Doutorado em Engenharia Química - Programa de Engenharia Química - Coppe - UFRJ.
 
 %Problema-exemplo Reator Semibatelada Schlegel(2004);


 % Opcoes do NLP
  optODE = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8 );
  optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-8),...
 'TolFun', 10^(-8), 'TolCon', 10^(-8), 'MaxFunEval', 10000);
  close all;
  
 % Horizonte de tempo....
 t0 = 0.0; % Tempo Inicial
 tf = 250; % Tempo Final
 x0 = [ 0.72; 0.05; 1]; %Seis estimativas iniciais
 wname = 'haar';  % tipo de base wavelets.
 tptr = 'minimaxi';
 sorh = 'h';
 scal = 'mln';
 is0 = 3;
 thrs = [];
for is = is0:10
    
   
% Numero de estagios: ns = 2, 4, 8, 16, and 32 ... 2^is
Cgrad=0;
grad=0;
Dgrad=0;
Pgrad=0;
if is==is0
   ns = (2^is);
    wL = 0*ones(ns,1);  % limite minimo da u
    wU = 0.001*ones(ns,1);  % limite maximo da u
    w0 = ones(ns,1)*0.0005;  % estimativa inicial da u...
    ts = [t0:(tf-t0)/ns:tf]; % Equiparticão padrão...
    point = [1:ns]*0;
end
 iu=is+2;
 % Chamando o FminCon...
 tic
 [wopt, erropt, iout] = fmincon( @(ws)obj(x0,ns,ts,ws,optODE), w0, [], [], [],...
 [], wL, wU, @(ws)ctr(x0,ns,ts,ws,optODE), optNLP);
 toc
 normGrad=norm(wopt(2:(length(wopt)-1)))/(length(wopt)-2);
 for ih=2:(length(wopt))
    Dgrad(ih-1)= abs(wopt(ih))-abs(wopt(ih-1));
 end
 il=1;
 
 for iw=1:length(Dgrad)
%      Dgsin=Dgrad(iw)*Dgrad(iw-1);
     if abs(Dgrad(iw))>0.2*normGrad
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
%  plot(tpq);
 
%  cA = appcoef(C,L,'db1',maxlev);
%  A = wrcoef('a',C,L,'db1',maxlev);
%  subplot(2,2,1); plot(A); 
%  title('Approximation')
%    for zl=1:maxlev
%     d = wrcoef('d',C,L,'db1',zl);
%     subplot(2,2,zl+1); plot(d); 
%     title('Detail')
%    end
   
 figure(iu+1000+530); % plotando os coeficientes...
 cfd = zeros(maxlev,length(wopt)); 
 for k = 1:maxlev
    d = detcoef(C,L,k); 
    d = d(ones(1,2^k),:); 
    cfd(k,:) = wkeep(d(:)',length(wopt)); 
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
 %Nc = norm(C);  % Norma dos coeficientes...
%  for i=1:(length(C))
%  if C(i)<0.000001;  %critério para eliminação de pontos...
%     C(i)=0.0;
%  end
%  end
%[thr,sorh,keepapp] = ddencmp('cmp','wv',wopt);  %thresholding
[w0,C,L] = wden(wopt,tptr,sorh,scal,maxlev,wname);
thr = thselect(wopt,tptr);
thrs(is)=thr;
%w0 = wdencmp('gbl',wopt,wname,maxlev,thr,sorh,keepapp);
 
[C,L] = wavedec(w0,maxlev,wname);
%tpq = wpdec(w0,maxlev,wname);
%  figure(iu+4000);
%plot(tpq);
figure(iu+2000+530);
w0 = waverec(C,L,wname); % transformada wavelets inversa...
cfd = zeros(maxlev,length(wopt)); 
 
for k = 1:maxlev  % Organizando o gráfico com os coeficientes de cada nivel
  d = detcoef(C,L,k); 
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
else
     w0=wopt;
     w0ns=wopt;
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
%             if (ind==0)
            w0ns(jj)=w0(j);      %nao vou dividir nesse caso...
            w0ns(jj+1)=w0(j+1); 
            ts(jj+1)=t0+((tf-t0)/(2*ns))*(2*j);  %
            ts(jj+2)=t0+((tf-t0)/(2*ns))*( (2*j+2) );  %
            point(jj)=1; point(jj+1)=1;
            jj=jj+2;
%            else
%               w0ns(jj)=w0(j);   w0ns(jj+1)=w0(j);    %aqui divido...
%               w0ns(jj+2)=w0(j+1); w0ns(jj+3)=w0(j+1);
%               ts(jj+1)=t0+((tf-t0)/(2*ns))*(2*(j-1)+1); ts(jj+2)=t0+((tf-t0)/(2*ns))*(2*(j-1)+2);
%               ts(jj+3)=t0+((tf-t0)/(2*ns))*(2*(j-1)+3); ts(jj+4)=t0+((tf-t0)/(2*ns))*(2*(j-1)+4);
%               point(jj)=0;point(jj+1)=0;point(jj+2)=0;point(jj+3)=0;
%               jj=jj+4; 
%            end
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
%                      point(jj)=0;point(jj+1)=0;point(jj+2)=0;point(jj+3)=0;
%                        ts(jj+1)=ts(jj)+(t0+((tf-t0)/(2.5*ns)));
%                        ts(jj)=ts(jj+1)-(t0+((tf-t0)/(2.5*ns)));
%                        ts(jj+2)=ts(jj+1)+(t0+((tf-t0)/(2.5*ns)));
%                        ts(jj+3)=ts(jj+2)+(t0+((tf-t0)/(2.5*ns)));
%                         ts(jj+4)=ts(jj+3)+(t0+((tf-t0)/(2.5*ns)));
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
 wL = 0*ones(ns,1);      % limite minimo da u
 wU = 0.001*ones(ns,1);  % limite maximo da u
 w0 = w0ns;
 %wL = -10*ones(ns,1);
 %wU = 10*ones(ns,1);
 end




 
  
  

 
 
 