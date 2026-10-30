
% This test illustrates that minimum MSE thresholds, outperform minimax,
% universal or FDR thresholds when we have the choice between soft- and hard
% thresholding
%
% In spite of ThreshLab's intelligent and fast routine for marginal 
% log-likelihood maximization in EBayesThresh approach, this routine is still
% the bottleneck in the simulation study when the zero inflated Laplace prior
% (with two parameters) is used. Simulations run much faster if only one
% parameter needs to be tested, so it is possible to make one parameter (a) to
% be "known". This speeds up the study, but also untruely favors EBayesThresh
% We can estimate a from the unobserved noise-free data and plug that value
% into ebayesthresh. Set afromnoisefree = 1; to run this option
%
% setup of Johnstone and Silverman, Ann.Stat., 2004, p. 1609
%
% alternatively, the simulation may adopt a true zero-inflated Laplacian prior 
% model for the generation of the noise-free components. Set truelaplace = 1;
% for this option

if exist('n') ~= 1, n = 10000; end
if exist('p') ~= 1, p = 0.005; end
K = p*n;
p = K/n;
if exist('mu0') ~= 1, mu0 = 7; end
if exist('nsimul') ~= 1, nsimul = 100; end
if exist('afromnoisefree') ~= 1, afromnoisefree = 0; end
if exist('truelaplace') ~= 1, truelaplace = 0; end

M = 500;
rand('state',M);
randn('state',M);

alfa = 0.05;

for simul= 1:nsimul,
   U = randn(n,1);

   if ~truelaplace,
      % next 2 lines are setup johnstone04:needles
      v = zeros(size(U));
      x = randperm(n); x = x(1:K); v(x) = mu0;
   else
      % next 3 lines are exact zero inflated Laplacian prior
      a = 1/mu0;
      r = rand(n,1);
      v = randlaplace(n,1,a).*(r<p);
   end

   W = v+U;


   stdevest = MAD1(W)/invcumgauss(3/4);
   vhatMSEsoft = minmsesoft(v,W);
   vhatSUREsoft = minsuresoft(W,stdevest);
   vhatGCVsoft = mingcvsoft(W);

   vhatMSEhard = minmsehard(v,W);
   vhatSUREhard = minsurehard(W,stdevest);
   vhatGCVhard = mingcvhard(W);


   if afromnoisefree,
      a = 1/mu0;
   else
      a = -1;
   end
   vhatEBthr = ebayesthresh(W,'prior','laplace','a',a,'thresrule','median');
   thrFDR = FDRthreshold(W,1,alfa);
   vhatFDR = HT(W,thrFDR);

   mseSUREsoft(simul) = norm(vhatSUREsoft - v)^2/n;
   mseGCVsoft(simul) = norm(vhatGCVsoft - v)^2/n;
   mseMSEsoft(simul) = norm(vhatMSEsoft - v)^2/n;

   mseSUREhard(simul) = norm(vhatSUREhard - v)^2/n;
   mseGCVhard(simul) = norm(vhatGCVhard - v)^2/n;
   mseMSEhard(simul) = norm(vhatMSEhard - v)^2/n;


   mseEBthr(simul) = norm(vhatEBthr - v)^2/n;
   mseFDR(simul) = norm(vhatFDR - v)^2/n;

simul

end

pos = (0:7)*0.6;
pos = (0:7);

verticalboxplot(mseSUREsoft,'linewidth',2,'color','b','pos',pos(1),'shaded')
hold on
verticalboxplot(mseGCVsoft,'linewidth',2,'color','r','pos',pos(2),'shaded')
verticalboxplot(mseMSEsoft,'linewidth',2,'color','g','pos',pos(3),'shaded')

verticalboxplot(mseFDR,'linewidth',2,'color','k','pos',pos(4))
verticalboxplot(mseEBthr,'linewidth',2,'color','m','pos',pos(5))

verticalboxplot(mseSUREhard,'linewidth',2,'color','b','pos',pos(6))
verticalboxplot(mseGCVhard,'linewidth',2,'color','r','pos',pos(7))
verticalboxplot(mseMSEhard,'linewidth',2,'color','g','pos',pos(8))

hold off

set(gca,'fontsize',14,'fontweight','b')

xtics = pos;
xlabs = {'SURE' 'GCV' 'MSE' 'FDR' 'EB' 'SURE' 'GCV' 'MSE'};
set(gca,'XTick',xtics)
set(gca,'XTickLabel',xlabs)
set(gcf,'position',[0 0 1000 500])
set(gcf, 'PaperPosition', [0 0 30 15]);

disp('Shaded boxplots are for SOFT thresholding, others are for hard')
disp('thresholding, or, in case of EB, for the EBayesThresh posterior median')
disp('threshold')
