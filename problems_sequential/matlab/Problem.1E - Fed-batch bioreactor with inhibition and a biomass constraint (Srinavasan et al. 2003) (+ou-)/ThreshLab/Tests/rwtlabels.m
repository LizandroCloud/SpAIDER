
% test procedure: image noise reduction by level dependent thresholding
% and non-decimated wavelet transform
% label data structure = given: determines which coefficients to remove

N = size(input,1);
M = size(label,1)-1;
% laatste rij labels komt overeen met schaalcoef.
wcoef = RWT2_BIO(input,M,ht,gt,sht,sgt);
wout = wcoef;

for j = 1:M
j
  % horizontal component

  W = zeros(size(input));
  X = W;
  W(:,:) = wcoef(j,1,:,:);
  X(:,:) = label(j,1,:,:);
  Wthr = W.*X;
  wout(j,1,:,:) = Wthr;

  % vertical subband

  W = zeros(size(input));
  X = W;
  W(:,:) = wcoef(j,2,:,:);
  X(:,:) = label(j,2,:,:);
  Wthr = W.*X;
  wout(j,2,:,:) = Wthr;

  % diagonal subband

  W = zeros(size(input));
  X = W;
  W(:,:) = wcoef(j,3,:,:);
  X(:,:) = label(j,3,:,:);
  Wthr = W.*X;
  wout(j,3,:,:) = Wthr;

end
alfa = 0.5;
output = IRT2_BIO(wout,h,g,sh,sg,alfa);
% noise2 = output - house;
% ratio2 = snr2(house,noise2);
image(output)

