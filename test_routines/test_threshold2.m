

% 'rigrsure' use principle of Stein's Unbiased Risk.
% 'heursure' is an heuristic variant of the first option.
% 'sqtwolog' for universal threshold sqrt(2*log(.)).
% 'minimaxi' for minimax thresholding.
 

 wname = 'db1';  % tipo de base wavelets.
 tptr = 'heursure';
 sorh = 'h';
 scal = 'sln';
 is0 = 2;
 thrs = [];
 
%  wgx=AA;
 
 
%   wopt=[1.22981903029542,1.33644862339731,1.43171055198538,1.51634322563408,1.59162603932803,1.65853028958594,1.71760791328755,1.76960136432644,1.8359388800008, 0.482386442306509,0.306097570742583,0.220828520500383,0.207826552552209,0.203979644457102,0.227571093042759,0.251103246197265,0.276487948584357,0,0,0.412676951057678,0.482888472175955,0,0.578319634121656,0.639321912762999 ];
%   wopt=[1.51124224601781,0,0,0,0,0,0,0.983623136013518,0,0.883761557525967,0,0.785578905934814,0.732179993265733,0.675441662755597,0.616183501762883,0.482386442306509,0.306097570742583,0.220828520500383,0.207826552552209,0.203979644457102,0.227571093042759,0.251103246197265,0.276487948584357,0,0,0.412676951057678,0.482888472175955,0,0.578319634121656,0.639321912762999,0.702347519033720,0.745418429101277,0.789326840053068,0.834079581125497,0.879422744893521,0.925166780431499,0.971444500295401,1.01827492708295,1.06540568104497,1.11265639458444,1.16017924479922,1.20780213349765,1.26953868355384,1.34586463364483,1.45382063098834,1.58443364180215,1.70192183440081,1.80686624101251,1.90062000682901,1.98427110221374,2.05843044697032,2.12400749090914,2.18230433708268,2.23438250522474,2.28095043755944,2.32283701309356,2.36108371296456,2.39655436674282,2.34172283782656,2.20496143083188,2.07207482108627,1.94293575885235,1.80239736890327,1.65388908164788,1.51833592982289,1.39408188752987,1.28046803130227,1.17658701239669,1.08101829560259,0.992757049058661,0.911432368576592,0.836474067204164,0.766893396825525,0.701980591376041,0.641452854859259,0.584887599206156,0.531517049079258,0.480763274682789,0.432371199954182,0.385973993016210,0.184946878867261];
%  wopt=[1.51124224601781,0,0,0,0,0,0,0.983623136013518,0,0.883761557525967,0,0.785578905934814,0,0.675441662755597,0.616183501762883,0,0.306097570742583,0.220828520500383,0.207826552552209,0,0,0.251103246197265,0.276487948584357,0,0,0.412676951057678,0.482888472175955,0,0.578319634121656,0,0,0.745418429101277,0.789326840053068,0.834079581125497,0,0,0.971444500295401,1.01827492708295,1.06540568104497,3,3,1.20780213349765,1.26953868355384,1.34586463364483,1.45382063098834,1.58443364180215,1.70192183440081,1.80686624101251,1.90062000682901,1.98427110221374,2.05843044697032,2.12400749090914,2.18230433708268,2.23438250522474,2.28095043755944,2.32283701309356,2.36108371296456,2.39655436674282,2.34172283782656,2.20496143083188,2.07207482108627,1.94293575885235,1.80239736890327,1.65388908164788,1.51833592982289,1.39408188752987,1.28046803130227,1.17658701239669,1.08101829560259,0.992757049058661,0.911432368576592,0.836474067204164,0.766893396825525,0.701980591376041,0.641452854859259,0.584887599206156,0.531517049079258,0.480763274682789,0.432371199954182,0.385973993016210,0.184946878867261]/8;
% wopt=wgx(1,:);
 wopt=[ 1 1 1 1 3 3 3 3 3 3 3 3 -1 -1 -1 -1 -1 -1 6 6 6 6 6 6 1 1 1 1 1 7 7 7];
x= randn(1,32); 
% wopt=[ 1  1 3 3 2 2 1 1 3 3 2 2 4 4 4 4 5 5 5 5 2 2 2 2 3 3 2 2 1 1 5 5 5 5 6 6 6 6 8 8 3 3 8 8 2 2 2 2 1 1 4 4 4 4 6 6 2 2 7 7 7 7 8 8];
% wopt = freqbrk;
sqrt_snr = 2;      % Set signal to noise ratio
init = 2055615866; % Set rand seed
% [xref,x] = wnoise(1,5,sqrt_snr,init);
% wopt=x;
woptx = wopt + x*0.8;
 maxlev = wmaxlev(length(woptx),wname);

[C,L] = wavedec(woptx,maxlev,wname);  % C sao os detalhes...

iesp = sum(C(L(1)+1:length(C)).^2);
esp = 1 + ( log(L(end)-L(1)) )^3/2 / (((L(end)-L(1)))^(1/2));
 tpq = wpdec(woptx,maxlev,wname);
 figure(43); % plotando os coeficientes...
 cfd = zeros(maxlev,length(woptx)); 
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
subplot(211), colormap(1-gray(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
title('Coeficientes da transformada wavelets') 
ylabel('Nivel')
j=1;
 Nc = norm(C(L(1)+1:length(C)))/length(C(L(1)+1:length(C)));% coefficient norm...
%  Nc = mean((C(L(2)+1:length(C))))
%  for i=1:(length(C))
%  if abs(C(i))<abs(Nc);  %critério para eliminação de pontos...
%     C(i)=0.0;
%  end
%  end
 [thr,sorh,keepapp] = ddencmp('cmp','wv',woptx);  %thresholding
 [Cth, y] = mingcvsoft(C);
%  Rvar=var(C);
%  sigmahat=median(C)/0.6745;
%  thr=bayes(C,sigmahat);
 Cth = sthresh(C,thr);
%  [thr, y] = mingcvhard(wopt);
% tt = (0:length(C)-1)/length(C)*thr;
% as= GCV1soft(C,tt)
%         [w0,C,L,s] = wden(wopt,tptr,sorh,scal,maxlev,wname);
%w0 = wdencmp('gbl',wopt,wname,maxlev,thr,sorh,keepapp);
 
% [C,L] = wavedec(w0,maxlev,wname);
%tpq = wpdec(w0,maxlev,wname);
%  figure(iu+4000);
%plot(tpq);
% figure(43);
d1 = detcoef(Cth,L,1);
w0 = waverec(Cth,L,wname); % transformada wavelets inversa...
cfd = zeros(maxlev,length(woptx)); 
chg = detcoef(Cth,L,2);
 Pgrad = Blips_test(0,w0,Cth,L)
for k = 1:maxlev  % Organizando o gráfico com os coeficientes de cada nivel
  d = detcoef(Cth,L,k); 
  d = d(ones(1,2^k),:); 
  cfd(k,:) = wkeep(d(:)',length(woptx)); 
end 

cfd = cfd(:); 
I = find(abs(cfd)<sqrt(eps)); 
cfd(I)=zeros(size(I)); 
cfd = reshape(cfd,maxlev,length(woptx));

% Plot discrete coefficients. 
subplot(212), colormap(1-gray(64)); 
img = image(flipud(wcodemat(cfd,64,'r'))); 
set(get(img,'parent'),'YtickLabel',[]); 
title('Coeficientes da transformada wavelets sem ruido') 
ylabel('Nivel')
 % dividindo estagios...
hold on;
figure (2)
%  stairs(w0)
hold on
stairs(wopt,'r--')
hold on;
stairs(woptx,'g')
