
function maxlines = connectmaxlines(maximage);

% connectmaxlines - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     maxlines = connectmaxlines(maximage);
%   Inputs
%     maximage  Binary image where each column corresponds to the location
%               of an observation and each row to a resolution level.
%               Pixel takes value 1 iff within this scale the WT reaches a
%               local maximum
%   Outputs
%     maxlines   Matrix where each column has the positions of a line of maxima
%                at successive scales. The number of columns equals the total
%                number of maxima lines.
%   Description
%      given a binary image indicating for each scale the locations of local
%      maxima, this procedure links those locations into lines of local maxima
%      across scales.
%   Background
%      The basic idea of this linking procedure is as follows:
%      Two maxima at successive scales are linked with each other if both are
%      the closest maximum to the other one. This procedure only links maxima
%      at immediately successive scales. In order to link maxima across gaps in
%      scales, use the procedure mergemaxlines.m
%   See also
%      help removelowmaxlines
%      help mergemaxlines
%      help deconnectmaxlines

M = size(maximage,1);

coarsemax = find(maximage(M,:));
C = length(coarsemax);
maxpos(M,1:C) = row(coarsemax);
coarsemaxlines = (1:C);
for j = M-1:-1:1,
   finemax = find(maximage(j,:));
   newmax = ones(size(finemax));
   for f = 1:length(finemax),
      [d closestcoarse(f)] = min(abs(finemax(f)-coarsemax)); end
   F = 0;
   for c = 1:length(coarsemax),
      [d f] = min(abs(coarsemax(c)-finemax)); 
      if closestcoarse(f) == c, % we have found a couple!
         F = F+1;
         % f = closestfine(c);
         newmax(f) = 0;
         l = coarsemaxlines(c);
         finemaxlines(f) = l;
         maxpos(j,l) = finemax(f);
      end
   end
   L = size(maxpos,2);
   for f = find(newmax),
      L = L+1;
      finemaxlines(f) = L;
      maxpos(j,L) = finemax(f);
   end
   coarsemaxlines = finemaxlines;
   clear finemaxlines;
   coarsemax = finemax;
end
hold off
for l = 1:size(maxpos,2),
   maxscales = find(maxpos(:,l));
   maxposl = maxpos(maxscales,l);
   plot(maxposl,maxscales)
   hold on
end
hold off
maxlines = maxpos;
