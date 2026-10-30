function maxlines = mergemaxlines(maxlines);

% mergemaxlines - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     maxlines = mergemaxlines(maxlines);
%   Inputs
%     maxlines   Matrix whose columns are plots of positions of local maxima as
%                function of scale
%   Outputs
%     maxlines   Modified version of input matrix. Dimensions may have changed.
%   Description
%      Lines of maxima across scales that act within overlapping intervals are
%      merged. This merging proceeds on a purely geometrical basis. The values
%      of the maxima are not involved.
%   See also
%      help connectmaxlines
%      help removelowmaxlines
%      help deconnectmaxlines

% This version uses relative distances in scales and positions.
% This implies that the result depends on the width of the interval we are
% working on. 

% Current version crashes on matrices with zero columns... could easily be
% dealt with by elimating those columns....wip


ll = size(maxlines,2);
smax = size(maxlines,1);
pmax = max(max(maxlines)) - min(min(maxlines(maxlines>0)))+1;
linedensity = sum(~~maxlines')/pmax;
reallines = ones(1,ll);
% We first find the minimum and maximum position/scale of each maxima line:
for l = 1:ll,
   maxscalesl = find(maxlines(:,l));
   maxxl = maxlines(maxscalesl,l);
   lineminpos(l) = min(maxxl);
   linemaxpos(l) = max(maxxl);
   lineminscale(l) = min(maxscalesl);
   linemaxscale(l) = max(maxscalesl);
   nlinepoints(l) = length(maxscalesl);
end
linemedianpos = (linemaxpos + lineminpos)/2;
% Next we check overlap among the x-values of different lines
overlap = checkoverlap(lineminpos,linemaxpos,lineminscale,linemaxscale);
% Now do the actual linking. Start with longest lines first. Note that a line
% from top to bottom will not be merged, as it cannot have any overlap.
linelength = linemaxscale - lineminscale + 1;
% [sortedlength orderoflines] = sort(-linelength);
[sortedlength orderoflines] = sort(-nlinepoints.*linemaxscale);
sortedlength = -sortedlength;
L = max(find(sortedlength>1));
for o = 1:L,
   l = orderoflines(o);
   candidates = find(overlap(l,:));
   % plotmaxlines(maxlines,'g')
   % hold on
   % plotmaxlines(maxlines(:,l),'r')
   % plotmaxlines(maxlines(:,candidates))
   % hold off
   % pause
   while ~isempty(candidates),
   % plotmaxlines(maxlines,'g')
   % hold on
   % plotmaxlines(maxlines(:,l),'r')
   % plotmaxlines(maxlines(:,candidates))
   % hold off
   % pause

      % find candidate with shortest jump
      jump = Inf;
      psps = zeros(length(candidates),4);
                 % psps = ( (positioo,scale) , (position,scale) )
                 % cooridinates of candidate links/bridges
      for c = 1:length(candidates), % compute begin and end point of bridge
         cc = candidates(c);
         if lineminscale(cc) < lineminscale(l),
            psps(c,2) = lineminscale(l);
            psps(c,4) = min(linemaxscale(cc),lineminscale(l));
         else
            psps(c,2) = linemaxscale(l);
            psps(c,4) = max(lineminscale(cc),linemaxscale(l));
         end
         psps(c,1) = maxlines(psps(c,2),l);
         psps(c,3) = maxlines(psps(c,4),cc);
      end
      % Now check if any of the candidate bridges could cross an existing line
      for c = 1:length(candidates),
         p1 = min(psps(c,1),psps(c,3));
         p2 = max(psps(c,1),psps(c,3));
         for s = min(psps(c,2),psps(c,4)):max(psps(c,2),psps(c,4)),
         if any(maxlines(s,:) < p2 & maxlines(s,:) > p1), psps(c,1) = Inf; end
         end
      end

      % jumps in scale directions are measured in a different way than jumps in
      % locations. Otherwise, there is a tendency to shoot beyond maxima on a
      % slightly different location:
      bridgelengths = abs(psps(:,4)-psps(:,2))/linelength(l) + ...
                      abs(psps(:,3)-psps(:,1)).*linedensity(psps(:,2))';
      % Do not allow jumps that are wider in scale than half of the length of
      % the existing line: (This is a very heuristical rule, so you might
      % consider to remove it, but it does prevent some really artificial
      % bridges) (added 10 feb 2005)
      c = abs(psps(:,4)-psps(:,2)) > linelength(l)/2 ;
      bridgelengths(c) = Inf;
      % in case two bridges have the same length, let the median position of
      % the involved lines of local maxima decide:
      bridgelengths = bridgelengths' + ...
                      abs(linemedianpos(l) - linemedianpos(candidates)).*...
                      linedensity(psps(:,2));
      [minlength c] = min(bridgelengths);
      if minlength < Inf,
         cc = candidates(c);
         s1 = lineminscale(l);
         s2 = linemaxscale(l);
         oldmaxline = maxlines(s1:s2,l);
         maxlines(1:s1,l) = maxlines(1:s1,cc);
         maxlines(s2:smax,l) = maxlines(s2:smax,cc);
         maxlines(s1:s2,l) = oldmaxline; 
                                      % this is to fill in the right values at
                                      % the end points of the original line.
         maxscalesl = find(maxlines(:,l));
         maxxl = maxlines(maxscalesl,l);
         lineminpos(l) = min(maxxl);
         linemaxpos(l) = max(maxxl);
         lineminscale(l) = min(maxscalesl);
         linemaxscale(l) = max(maxscalesl);
         maxlines(:,cc) = 0;
         lineminpos(cc) = 0;
         linemaxpos(cc) = 0;
         lineminscale(cc) = 0;
         linemaxscale(cc) = 0;
         reallines(cc) = 0;
         ccc = [candidates find(overlap(cc,:))];
         overlap(ccc,ccc) =  checkoverlap(lineminpos(ccc),linemaxpos(ccc),...
                                          lineminscale(ccc),linemaxscale(ccc));
         overlap(:,cc) = 0;
         overlap(cc,:) = 0;
         candidates = ccc(find(overlap(l,ccc)));
      else, candidates = []; end
   end
end
maxlines = maxlines(:,find(reallines));


function overlap = checkoverlap(lineminpos,linemaxpos,lineminscale,linemaxscale);

% This function finds close lines of maxima that are candidates to be linked
% into a single line. Generally speaking, if there is overlap in locations,
% linking is promoted. If there is overlap in scales, linking is discouraged.

% first, add some additional uncertainty bands (It may happen that two lines 
% do not actually overlap, but are nevertheless to be linked)
lineminpos1 = 2*lineminpos - linemaxpos+1;
linemaxpos = 2*linemaxpos - lineminpos-1;
lineminpos = lineminpos1;

ll = length(lineminpos);
overlap = zeros(ll,ll);
overlapscale = zeros(ll,ll);
for l = 2:ll, 
   % first look at range of positions
   nooverlappos = lineminpos(l) > linemaxpos(1:l-1) | ...
                  linemaxpos(l) < lineminpos(1:l-1);
   % next look at range of scales: if the two lines have an overlap in scale
   % which is larger than the non-overlapping parts, there is no need to jump
   overlapscalel = min([linemaxscale(l)-lineminscale(1:l-1) ; ...
                        linemaxscale(1:l-1)-lineminscale(l)]);
   % This value is negative for non-overlapping lines. It overestimates the
   % overlap for a line that fully overlaps another, but this overestimation
   % will only reinforce the decision not to make these lines candidates for
   % mutual linking (the higher the value of overlap, the more unlikely that
   % the lines are going to be linked)
   nonoverlapscalel = linemaxscale(l) - lineminscale(l) - overlapscalel;
   nonoverlapscalem = linemaxscale(1:l-1) - lineminscale(1:l-1) - overlapscalel;
   dominatedoverlapscalel = (overlapscalel > nonoverlapscalel/4) | ...
                           (overlapscalel > nonoverlapscalem/4);
   m = find(~nooverlappos & ~dominatedoverlapscalel); 
   overlap(l,m) = 1;
   overlap(m,l) = 1;
end
