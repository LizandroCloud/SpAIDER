

% 'rigrsure' use principle of Stein's Unbiased Risk.
% 'heursure' is an heuristic variant of the first option.
% 'sqtwolog' for universal threshold sqrt(2*log(.)).
% 'minimaxi' for minimax thresholding.
 

 wname = 'db1';  % tipo de base wavelets.
 tptr = 'sqtwolog';
 sorh = 'h';
 scal = 'sln';
 is0 = 2;
 thrs = [];
wopt=[ 1 1 1 1 3 3 3 3 3 3 3 3 -1 -1 -1 -1 -1 -1 6 6 6 6 6 6 1 1 1 1 1 7 7 7];

%  wopt=[1.0000E+00
% 1.0000E+00
% 1.0000E+00
% 9.2248E-01
% 8.4496E-01
% 7.9845E-01
% 7.3643E-01
% 6.8217E-01
% 6.3566E-01
% 5.8915E-01
% 5.5039E-01
% 5.1163E-01
% 4.7287E-01
% 4.4186E-01
% 4.1085E-01
% 3.7984E-01
% 3.4109E-01
% 3.1008E-01
% 2.7907E-01
% 2.5581E-01
% 2.3256E-01
% 2.1705E-01
% 1.9380E-01
% 1.7829E-01
% 1.7054E-01
% 1.6279E-01
% 1.5504E-01
% 6.5116E-01
% 6.5116E-01
% 6.4341E-01
% 6.3566E-01
% 6.2791E-01
% ]';
% 
% wopt = [7.6154E-01
% 7.0769E-01
% 6.6154E-01
% 6.2308E-01
% 5.8462E-01
% 5.5385E-01
% 5.2308E-01
% 4.9231E-01
% 4.6923E-01
% 4.4615E-01
% 4.3077E-01
% 4.0769E-01
% 3.8462E-01
% 3.7692E-01
% 3.6154E-01
% 3.3846E-01
% 3.2308E-01
% 3.0769E-01
% 2.8462E-01
% 2.6923E-01
% 2.6923E-01
% 2.5385E-01
% 2.5385E-01
% 2.4615E-01
% 2.4615E-01
% 2.4615E-01
% 2.4615E-01
% 2.4615E-01
% 2.4615E-01
% 2.4615E-01
% 2.3846E-01
% 2.3846E-01]';

%     x = wgn(1,32,0);
for ki=1:10
woptx = wopt + x*ki*0.1
 maxlev = wmaxlev(length(woptx),wname);
[C,L] = wavedec(woptx,maxlev,wname);  % C sao os detalhes...
 tpq = wpdec(woptx,maxlev,wname);
 figure(55); % plotando os coeficientes...
 cfd = zeros(maxlev,length(wopt)); 
 for k = 1:maxlev
    d = detcoef(C,L,k); 
    d = d(ones(1,2^k),:); 
    cfd(k,:) = wkeep(d(:)',length(woptx)); 
 end 

 cfd = cfd(:); 
 I = find(abs(cfd)<sqrt(eps)); 
 cfd(I)=zeros(size(I)); 
 cfd = reshape(cfd,maxlev,length(woptx));

% Plotando os coeficientes. 
subplot(212), colormap(1-pink(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
% title('Coeficientes da transformada wavelets') 
ylabel('Wavelet level')
j=1;
% trr = norm(a)/length(a);
%  Nc = norm(C);  % Norma dos coeficientes...
%  for i=1:(length(C))
%  if C(i)<trr;  %critério para eliminação de pontos...
%     C(i)=0.0;
%  end
%  end
 
 
 
%[thr,sorh,keepapp] = ddencmp('cmp','wv',wopt);  %thresholding
[w0,Cth,Lth,s,thr] = spaider_wden(woptx,tptr,sorh,scal,maxlev,wname);
% thr = thselect(woptx,tptr);
% w0 = wdencmp('gbl',woptx,wname,maxlev,thr,sorh,keepapp);
%  
%  [C,L] = wavedec(w0,maxlev,wname);
tpq = wpdec(w0,maxlev,wname);
%   figure(iu+4000);
% plot(tpq);
figure(25);
w0 = waverec(Cth,L,wname); % transformada wavelets inversa...
cfd = zeros(maxlev,length(woptx)); 
 
for k = 1:maxlev  % Organizando o gráfico com os coeficientes de cada nivel
  d = detcoef(Cth,L,k); 
  d = d(ones(1,2^k),:); 
  cfd(k,:) = wkeep(d(:)',length(w0)); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(w0));

% Plot discrete coefficients. 
subplot(211), colormap(1-pink(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
% title('Coeficientes da transformada wavelets sem ruido') 
ylabel('Wavelet level')
 % dividindo estagios...
hold on;
figure (ki)
%  stairs(w0)
hold on
subplot(212)
stairs(wopt,'r')
hold on;
stairs(woptx,'b')
hold on;
stairs(w0,'k-')

mse(ki) = norm(wopt-w0)
ylabel('Dataset values')
xlabel('Sampled data')

legend('Clean dataset','Full dataset','Sure')



peaksnr(ki) = psnr(w0,wopt);
end
hold on;
figure(99)
subplot(311)
hold on;
plot(mse,'o')
xlabel('s')
ylabel('MSE')
subplot(312)
plot(x,'k')

hold on;
subplot(313)
plot(peaksnr,'x')
xlabel('s')
ylabel('MSE')
