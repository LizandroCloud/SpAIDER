
function optimalunbalanced = findallunbalancedhaar(w,M,maxlines,y,pdf)
         
% optimalunbalanced  - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     optimalunbalanced = findallunbalancedhaar(w,M,maxlines,y,pdf);
%   Inputs
%     w         Continuous Haar or Haar-Fisz wavelet coefficients
%     M         number of dyadic scales
%     maxlines  selected lines of local maxima across scales in w
%     y         observations
%     pdf       string: 'normal' or 'poisson' for Haar or Haar-Fisz analysis
%             respectively
%   Outputs
%     optimalunbalanced =
%     {optimalscale,optimalposition,optimalleft,optimalright,optimalcoef}
%     all of them vectors of length equal to number of maxlines
%     for every line of local maxima, the output includes
%         optimalscale        the most significant scale
%         optimalposition     the most significant position
%         optimalleft         the most significant left scale
%         optimalright        the most significant right scale
%         optimalcoef         the value of significance
%       
%    See also
%      help findoptimalunbalancedhaar: routine called for every maxline

if nargin < 5, pdf = 'normal'; end

% We find the scale with maximum value on each line. 

% we need on each scale the number of observations involved in the computation
% of the coefficients. We use the following trick

Mm = size(w,1); m = Mm-M;
[dummy LR] = CHT(1,M,m);
LR = LR/2;
% LR contains for every scale j the number of observations on the left hand
% side and the number of observations on the right hand side of the balanced
% Haar wave.
for line = 1:size(maxlines,2),
   maxscales = find(maxlines(:,line));
   maxx = maxlines(maxscales,line);
   clear maxw
   for i = 1:length(maxscales), maxw(i) = w(maxscales(i),maxx(i)); end
   [wi i] = max(abs(maxw));
   k = maxx(i);
   j = maxscales(i);
   scale = LR(j);
   kmin = min(maxx);
   kmax = max(maxx);

   [k,l,r,wu] = findoptimalunbalancedhaar(y,k,kmin,kmax,scale,'poisson');

   optimalscale(line) = j;
   optimalleft(line) = l;
   optimalright(line) = r;
   optimalcoef(line) = wu;
   optimalposition(line) = k;
end
optimalunbalanced = {optimalscale,optimalposition,optimalleft,optimalright,...
                     optimalcoef};
