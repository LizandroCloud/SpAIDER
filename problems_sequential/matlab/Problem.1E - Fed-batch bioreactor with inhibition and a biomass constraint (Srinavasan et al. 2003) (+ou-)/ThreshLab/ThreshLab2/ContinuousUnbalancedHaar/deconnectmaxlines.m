function maxlines = deconnectmaxlines(maxlines);

% deconnectmaxlines - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     maxlines = deconnectmaxlines(maxlines);
%   Inputs
%     maxlines   Matrix whose columns are plots of positions of local maxima as
%                function of scale
%   Outputs
%     maxlines   Modified version of input matrix. Dimensions may have changed.
%   Description
%      This routine is the third in a sequence of four routines that
%      constructs lines of local maxima across scales:
%        1) connectmaxlines
%        2) removelowmaxlines
%        3) deconnectmaxlines
%        4) mergemaxlines
%      The first routine actually constructs lines of maxima from a binary
%      image
%      The second routine eliminates lines of local maxima that never reach
%      signficance
%      This, third routine scissors a single line into two separate lines if a
%      jump in location between two adjacent scales is too large
%      The fourth routine joins lines that are very/too close into a single
%      line
%   See also
%      help connectmaxlines
%      help mergemaxlines
%      help removelowmaxlines

ll = size(maxlines,2);
L = ll;
smax = size(maxlines,1);
reallines = ones(1,ll);
l = 1;
for l = 1:L,
% while l <= ll,
   maxscalesl = find(maxlines(:,l));
   maxxl = maxlines(maxscalesl,l);
   r = length(maxxl);
   if r > 1,
      jumps = (maxxl(2:r) - maxxl(1:r-1) ) ./ ...
              (maxscalesl(2:r) - maxscalesl(1:r-1) );
   else
      jumps = 0;
   end
   stdevjumps = sqrt(var2(jumps));
   % stdevjumps = MAD1(jumps)/invcumgauss(3/4);
   B = min(20,floor(3*stdevjumps));
   while ismember(B,abs(jumps)), B = B+1; end
   meanjumps = mean(jumps);
   if length(jumps) == 1, meanjumps = 0; B = 20; end
   outliers = find(abs(jumps-meanjumps) > B);
   outliers = sort(outliers);
   for o = 1:length(outliers),
      oo = outliers(o);
      s = maxscalesl(oo);
      ll = ll+1;
      maxlines = [maxlines zeros(smax,1)];
      maxlines(1:s,ll) = maxlines(1:s,l);
      maxlines(1:s,l) = 0;
   end
   l = l+1;
end
