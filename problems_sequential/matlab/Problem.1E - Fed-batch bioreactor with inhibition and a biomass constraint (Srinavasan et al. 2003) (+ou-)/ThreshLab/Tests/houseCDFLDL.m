
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
 noise1 = whitenoise;
% or
% noise1 = filter2(H,whitenoise);
% following lines for proper rescaling of the noise, such that desired
% signal-to-noise level is obtained.
noise1 = double(noise1);
noise1 = noise1 - mean2(noise1);
aim = 3;
ratio1 = snr2(house,noise1);
factor = 10^((ratio1 - aim)/20);
noise1= noise1*factor;
ratio1 = snr2(house,noise1);

input = house+noise1;
image(input);
colormap(housemap);
noisefree = house;

[h ht] = CDF_LDL(4,4);                               
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht)

M = 3;
alfa = 0.5;
