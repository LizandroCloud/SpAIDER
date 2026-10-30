

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
wopt=a;
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
trr = norm(a)/length(a);
 Nc = norm(C);  % Norma dos coeficientes...
 for i=1:(length(C))
 if C(i)<trr;  %critério para eliminação de pontos...
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

