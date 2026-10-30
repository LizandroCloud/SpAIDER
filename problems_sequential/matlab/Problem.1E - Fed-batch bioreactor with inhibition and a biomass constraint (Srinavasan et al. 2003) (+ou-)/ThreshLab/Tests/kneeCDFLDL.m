
knee = imread('/home/maarten/public_html/images/knee.tif','tif');
knee = double(knee);
kneemap = [(0:255)' (0:255)' (0:255)']/255;
image(knee)
colormap(kneemap)

input = knee;

[h ht] = CDF_LDL(4,4);
[h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);

M = 2;
alfa = 0.5;
