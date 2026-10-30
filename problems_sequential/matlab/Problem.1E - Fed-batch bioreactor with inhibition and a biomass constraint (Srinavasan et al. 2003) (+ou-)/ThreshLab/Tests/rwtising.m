
% test procedure: image noise reduction by level dependent thresholding
% and non-decimated wavelet transform (+ GCV - estimation)
% improved by MRF-model 
% WARNING: this may take a while of computation!

N = size(input,1);
wcoef = RWT2_BIO(input,M,ht,gt,sht,sgt);
wout = wcoef;
label = zeros(size(wcoef));
maxsample = 10;

for j = 1:M
j
  % horizontal component

  W = zeros(size(input));
  W(:,:) = wcoef(j,1,:,:);
  [Wthr, thr] = mingcv2(W);
  X = (abs(W) > thr);
  NX = size(X,1)*size(X,2);
  sigsig = var2(Wthr-W);
  sig = sqrt(sigsig);
  Wmax = max(max(abs(W)));
  mu = ((NX + 1)*Wmax - sig)/NX;
muvec(j,1) = mu;
sigvec(j,1) = sig;
  X0 = medfilt2(X,[3 3]);
  tau = MPLEising(X0);
  % tau = 0.5;
  Y = X;
  disp('Metropolis sampler...')
  for sample = 1:maxsample
     X = metropolisIsing(X,W,sig,tau,thr);
     Y = Y + X;
     sample
  end
  X = Y/(maxsample+1);
  Wthr = W.*X;
  wout(j,1,:,:) = Wthr;
  label(j,1,:,:) = X;
  gcvthr(j,1) = thr;
  mrftau(j,1) = tau;

  % vertical subband

  W = zeros(size(input));
  W(:,:) = wcoef(j,2,:,:);
  [Wthr, thr] = mingcv2(W);
  X = (abs(W) > thr);
  NX = size(X,1)*size(X,2);
  sigsig = var2(Wthr-W);
  sig = sqrt(sigsig);
  Wmax = max(max(abs(W)));
  mu = ((NX + 1)*Wmax - sig)/NX;
muvec(j,2) = mu;
sigvec(j,2) = sig;
  X0 = medfilt2(X,[3 3]);
  tau = MPLEising(X0);
  % tau = 0.5;
  Y = X;
  disp('Metropolis sampler...')
  for sample = 1:maxsample
     X = metropolisIsing(X,W,sig,tau,thr);
     Y = Y + X;
     sample
  end
  X = Y/(maxsample+1);
  Wthr = W.*X;
  wout(j,2,:,:) = Wthr;
  label(j,2,:,:) = X;
  gcvthr(j,2) = thr;
  mrftau(j,2) = tau;

  % diagonal subband

  W = zeros(size(input));
  W(:,:) = wcoef(j,3,:,:);
  [Wthr, thr] = mingcv2(W);
  X = (abs(W) > thr);
  NX = size(X,1)*size(X,2);
  sigsig = var2(Wthr-W);
  sig = sqrt(sigsig);
  Wmax = max(max(abs(W)));
  mu = ((NX + 1)*Wmax - sig)/NX;
muvec(j,3) = mu;
sigvec(j,3) = sig;
  X0 = medfilt2(X,[3 3]);
  tau = MPLEising(X0);
  % tau = 0.5;
  Y = X;
  disp('Metropolis sampler...')
  for sample = 1:maxsample
     X = metropolisIsing(X,W,sig,tau,thr);
     Y = Y + X;
     sample
  end
  X = Y/(maxsample+1);
  Wthr = W.*X;
  wout(j,3,:,:) = Wthr;
  label(j,3,:,:) = X;
  gcvthr(j,3) = thr;
  mrftau(j,3) = tau;

end
alfa = 0.5;
output = IRT2_BIO(wout,h,g,sh,sg,alfa);
% noise2 = output - house;
% ratio2 = snr2(house,noise2);
image(output)

