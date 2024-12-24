% euler1   Programa para o calculo da E.D.O. y' = x + y 
%          Metodo de Euler 
%          Condicao de contorno: y(0) = 1 
clear; 
% Condicao de contorno 
x(1) = 0;    
y(1) = 1;    
n = input('Numero de intervalos: '); 
xf = input('Valor de x final: '); 
h = (xf - x(1))/n; 
for i = 1:n 
    f(i) = x(i) + y(i); 
    x(i+1) = x(i) + h; 
    y(i+1) = y(i) + h*f(i); 
end 
% Calculo da solucao exata 
xe = 0:(xf/n):xf; 
ye = 2*exp(xe) - xe - 1; 
e = norm(ye-y)
% Grafico comparativo: solucao pelo metodo de Euler e a solucao exata 
plot(xe,ye,'-r',x,y,'ob');  
xlabel('x'); ylabel('y = y(x)'); 
% legend('Exato','Euler'); 

% 'rigrsure' use principle of Stein's Unbiased Risk.
% 'heursure' is an heuristic variant of the first option.
% 'sqtwolog' for universal threshold sqrt(2*log(.)).
% 'minimaxi' for minimax thresholding.
 

 wname = 'haar';  % tipo de base wavelets.
 tptr = 'heursure';
 sorh = 'h';
 scal = 'mln';
 is0 = 2;
 thrs = [];
 wopt=y;
 maxlev = WMAXLEV(length(wopt),wname);
[C,L] = wavedec(wopt,maxlev,wname);  % C sao os detalhes...
 tpq = wpdec(wopt,maxlev,wname);
 figure(55); % plotando os coeficientes...
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
trr = norm(y)/length(y);
 Nc = norm(C);  % Norma dos coeficientes...
 for i=1:(length(C))
 if norm(C(i))<trr;  %critério para eliminação de pontos...
    C(i)=0.0;
 end
 end
% %[thr,sorh,keepapp] = ddencmp('cmp','wv',wopt);  %thresholding
% [w0,C,L] = wden(wopt,tptr,sorh,scal,maxlev,wname);
% thr = thselect(wopt,tptr);
% %w0 = wdencmp('gbl',wopt,wname,maxlev,thr,sorh,keepapp);
%  
% [C,L] = wavedec(w0,maxlev,wname);
%tpq = wpdec(w0,maxlev,wname);
%  figure(iu+4000);
%plot(tpq);
figure(25);
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
 figure(1000);
 plot(w0,'o');hold on; plot(y);
 jj=1;
 t0=0;
 ns = length(w0);
 ts = w0/length(w0) + 1;
 tf = length(w0);
 
 j =1;
 for i=1:1:ns-1
     if w0(i)==w0(i+1)
         w00(j)=w0(i+1);
         j = j+1;
     else
         w00(j)=w0(i);
     end
 end
         
         