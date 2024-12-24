
% This test illustrates the use of GCV (with several alternative definitions)
% in selecting the optimal threshold in a hard threshold variable selection
% scheme.
%
% observational model: Y = beta + noise
% 
% with beta from zero-inflated Laplace distribution with P(beta_i=0) = 1-p
% p is degreeofsparsity.
%
% the following parameters can be specified (with default values)
%
%     realrandom: (default 0) binary value; if 0, then pseudo random numbers
%                  are always the same.
%     degreeofsparsity (default 0.05) value p defined above
%     samplesize (default 2000)

if exist('realrandom') ~= 1, realrandom = 0; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if exist('samplesize') ~= 1, samplesize = 2000; end
% 0.35

if ~realrandom,
   randn('state',2000);
   rand('state',2000);
end

n = samplesize;
nthr = 100;
p = degreeofsparsity;

a = 1/5;
r = rand(1,n);
beta = randlaplace(1,n,a).*(r<p);

noise = randn(size(beta));

Y = beta+noise;

% Y = 2*Y; noise = 2*noise; beta = 2*beta;

stdev = sqrt(var2(noise));
stdevestMAD = MAD1(Y)/invcumgauss(0.75);
stdevestGCV = sqrt(varestGCV(Y));

tuniv = sqrt(2*log(n))*stdev;

tmax = tuniv;

tt = (0:nthr)/nthr*tmax;

figure(1)
plot(tt,GCV1soft(Y,tt),'r','linewidth',2)
aa = axis; aa(3) = 0;
plot(tt,MSE1thresh(beta,Y,tt,1),'linewidth',2)
hold on
plot(tt,SURE1soft(Y,tt,stdev),'m','linewidth',2)
plot(tt,GCV1soft(Y,tt),'r','linewidth',2)
plot(tt,GCV1soft(Y,tt)-stdev^2,'r--','linewidth',2)
hold off
axis(aa)
title('Soft thresholding','fontsize',14,'fontweight','bold')
legend('MSE','SURE','GCV','GCV-\sigma^2')
set(gca,'fontsize',14,'fontweight','bold')

[Ythr thr] = mingcvsoft(Y);

tt = 2*tt;
figure(2)
plot(tt,MSE1thresh(beta,Y,tt,0),'linewidth',2)
hold on
plot(tt,GCV1hard(Y,tt,stdevestMAD),'r','linewidth',2)
plot(tt,GCV1hard(Y,tt,stdevestGCV),'k','linewidth',2)
plot(tt,GCV1hard(Y,tt,stdev),'m','linewidth',2)
plot(tt,GCV1softhard(Y,tt),'g','linewidth',2)
hold off
title('Hard thresholding; GCV via soft thr.','fontsize',14,'fontweight','bold')
legend('MSE','GCV with \sigma by MAD/\Phi^{-1}(0.75)',...
       'GCV with \sigma by GCV','GCV with exact \sigma','naive GCV')
set(gca,'fontsize',14,'fontweight','bold')
aa = axis;
aa(2) = 2*tmax;
axis(aa)

figure(3)
plot(tt,MSE1thresh(beta,Y,tt,0),'linewidth',2)
hold on
plot(tt,GCV1hardDoFmirror(Y,tt,stdevestMAD),'r','linewidth',2)
plot(tt,GCV1hardDoFmirror(Y,tt,stdevestGCV),'k','linewidth',2)
plot(tt,GCV1hardDoFmirror(Y,tt,stdev),'m','linewidth',2)
hold off
title('Hard thresholding; GCV with corrected d.o.f.',...
       'fontsize',14,'fontweight','bold')
legend('MSE','GCV with \sigma by MAD/\Phi^{-1}(0.75)',...
       'GCV with \sigma by GCV','GCV with exact \sigma')
set(gca,'fontsize',14,'fontweight','bold')
axis(aa)

figure(4)
plot(tt,MSE1thresh(beta,Y,tt,0),'linewidth',2)
hold on
plot(tt,SURE1hard(Y,tt,stdev),'m','linewidth',2)
plot(tt,GCV1hard(Y,tt,stdev)-stdev^2,'k','linewidth',2)
plot(tt,GCV1hardDoFmirror(Y,tt,stdev)-stdev^2,'r','linewidth',2)
hold off
title('GCV versus SURE','fontsize',14,'fontweight','bold')
legend('MSE','SURE','GCV-\sigma^2 via ST','GCV-\sigma^2; mirror corr. d.o.f.')
set(gca,'fontsize',14,'fontweight','bold')
axis(aa)





[Ythr opthrHT] = mingcvhard(Y);
opthrHT
[Ythr opthrHT] = mingcvhard(Y,'maxthr',max(tt));
opthrHT
