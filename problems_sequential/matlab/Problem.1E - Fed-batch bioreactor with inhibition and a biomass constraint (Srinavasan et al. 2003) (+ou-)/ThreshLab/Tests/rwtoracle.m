
% test procedure: image noise reduction by level dependent thresholding
% and non-decimated wavelet transform (oracle thresholds)

N = size(input,1);
wcoef = RWT2_BIO(input,M,ht,gt,sht,sgt);
vcoef = RWT2_BIO(noisefree,M,ht,gt,sht,sgt);
wout = wcoef;
label = zeros(size(wcoef));

for j = 1:M
j
  % horizontal component

  w = zeros(size(input));
  w(:,:) = wcoef(j,1,:,:);
  v = zeros(size(input));
  v(:,:) = vcoef(j,1,:,:);
  [wt, thr] = oracle(v,w);
  wout(j,1,:,:) = wt;
  label(j,1,:,:) = (abs(v) > thr);
  oraclethr(j,1) = thr;

  % diagonal subband

  w = zeros(size(input));
  w(:,:) = wcoef(j,2,:,:);
  v = zeros(size(input));
  v(:,:) = vcoef(j,2,:,:);
  [wt, thr] = oracle(v,w);
  wout(j,2,:,:) = wt;
  label(j,2,:,:) = (abs(v) > thr);
  oraclethr(j,2) = thr;

  % vertical subband

  w = zeros(size(input));
  w(:,:) = wcoef(j,3,:,:);
  v = zeros(size(input));
  v(:,:) = vcoef(j,3,:,:);
  [wt, thr] = oracle(v,w);
  wout(j,3,:,:) = wt;
  label(j,3,:,:) = (abs(v) > thr);
  oraclethr(j,3) = thr;

end
alfa = 0.5;
output = IRT2_BIO(wout,h,g,sh,sg,alfa);
image(output)

