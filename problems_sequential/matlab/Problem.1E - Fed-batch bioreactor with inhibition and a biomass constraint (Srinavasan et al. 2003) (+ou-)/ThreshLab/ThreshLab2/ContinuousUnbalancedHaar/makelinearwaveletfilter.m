
function hhthh = makelinearwaveletfilter(type,n,nt,M)

% makelinearwaveletfilter - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     hhthh = makelinearwaveletfilter(type,n,nt,M);
%   Inputs
%     type       (ort)  'Haar','Coiflet','Daubechies','Symmlet'
%     type       (bi-ort) 'CDF', 'CDFless' (CDF with less dissimilar lengths)
%     n          number of (primal)  vanishing moments
%     nt         number of dual vanishing moments
%     M          number of resolution levels
%   Outputs
%     hhthh      linear filter
%   Description
%     This procedure constructs a linear, time-invariant filter corresponding
%     to a translation-invariant wavelet analysis followed by a reconstruction
%     from its coarse scaling coefficients only, i.e., the following algorithm:
%     
%     [h ht] = makeWTfilter(type,n,nt);
%     [h, g, ht, gt, sh, sg, sht, sgt] = mirrorshifts(h,ht);
%     Wt = RWT_BIO(eye(N),M,ht,gt,sht,sgt);
%         (this is not possible as such, since RWT_BIO does not run on
%         matrices)
%     XWt = zeros(size(Wt)); XWt(scaling coeff) = Wt(scaling coeff);
%     WXWt = IRT_BIO(XWt,h,g,sh,sg);

if nargin == 3, M = nt; nt = n; end
if nargin == 2, M = n; nt = n; end
[h ht] = makeWTfilter(type,n,nt);

ht0 = ht(length(ht):-1:1);
hht = ht0;
for j = 2:M,
   ht00 = zeros(1,length(ht0)*2-1);
   ht00(1:2:length(ht00)) = ht0;
   hht = conv(ht00,hht);
   ht0 = ht00;
end

h0 = h;
hh = h;
for j = 2:M,
   h00 = zeros(1,length(h0)*2-1);
   h00(1:2:length(h00)) = h0;
   hh = conv(h00,hh);
   h0 = h00;
end

hhthh = conv(hht,hh);
hhthh = hhthh/2^M;
