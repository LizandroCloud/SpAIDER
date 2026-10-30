
% test procedure: image noise reduction by level dependent thresholding
% and non-decimated wavelet transform (+ GCV - estimation)

N = size(input,1);
wcoef = RWT2_BIO(input,M,ht,gt,sht,sgt);
wout = wcoef;
label = zeros(size(wcoef));
clear gcvthr

for j = 1:M
j
  % horizontal component

  x = zeros(size(input));
  x(:,:) = wcoef(j,1,:,:);
  [y, thr] = mingcv2(x);
  wout(j,1,:,:) = y;
  label(j,1,:,:) = (abs(x) > thr);
  gcvthr(j,1) = thr;

  % vertical subband

  x = zeros(size(input));
  x(:,:) = wcoef(j,2,:,:);
  [y, thr] = mingcv2(x);
  wout(j,2,:,:) = y;
  label(j,2,:,:) = (abs(x) > thr);
  gcvthr(j,2) = thr;

  % diagonal subband

  x = zeros(size(input));
  x(:,:) = wcoef(j,3,:,:);
  [y, thr] = mingcv2(x);
  wout(j,3,:,:) = y;
  label(j,3,:,:) = (abs(x) > thr);
  gcvthr(j,3) = thr;

end
alfa = 0.5;
output = IRT2_BIO(wout,h,g,sh,sg,alfa);
image(output)

