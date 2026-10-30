
house = imread('~/Matlab/ThreshLab/Tests/TestImages/house.tif','tif');
house = double(house);
housemap = [(0:255)' (0:255)' (0:255)']/255;
image(house)
colormap(housemap)
whitenoise = randn(size(house));

% make 2D highpass-filter for noise coloring
% fir1: 1D filter-design (signal processing toolbox)
% help fir1 for more info
B = fir1(100,0.6,'high');
H = ftrans2(B); % from image processing toolbox
% chose colored or white noise:
% noise1 = whitenoise;
% or
% noise1 = filter2(H,whitenoise);
% following lines for proper rescaling of the noise, such that desired
% signal-to-noise level is obtained.
aim = 10;
ratio1 = snr2(house,noise1);
factor = 10^((ratio1 - aim)/20);
noise1= noise1*factor;
ratio1 = snr2(house,noise1);

input = house+noise1;
image(input);
colormap(housemap);
noisefree = house;

% qmf = MakeONFilter('Daubechies',4)
% h = qmf;
% ht = h;
% clear g
% p = length(h);
% for i = 1:p
%   g(i) = h(p-i+1)*(-1)^(i+1);
% end
% gt = g;
% sh = 0;
% sg = 2-p;
% sht = 0;
% sgt = 2-p;

h = Haar;
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h);
M = 3;
alfa = 0.5;
