
% This routine is called by the routine searchformultiplechangepoints.m
%
% Its starts from the analysis on all selected lines of local maxima,
% summarised in
% optimalunbalanced =
% {optimalscale,optimalposition,optimalleft,optimalright, optimalcoef};
%
% This routine makes the final selection of change points.
%
% It does so in three different ways, and the user may choose among those three
% possibilities
%
% 1) cpaic are change points selected by minimizing AIC. ffaic is the
%    corresponding piecewise constant fit of the data y
% 2) cpdaic is the largest set of change points selected such that adding one
%    after the other always causes a decrease of AIC (so, this is selection of
%    the smallest local minimum of the AIC curve. ffaic is the corresponding
%     piecewise constant fit of the data y
% 3) cpoutliersdaic is a set of change points that cause a significant decrease
%    in AIC value. All other decreases in AIC are considered insignificant and
%    the corresponding change points are not taken into account.

optimalcoef     = optimalunbalanced{5};
optimalright    = optimalunbalanced{4};
optimalleft     = optimalunbalanced{3};
optimalposition = optimalunbalanced{2};
optimalscale    = optimalunbalanced{1};
aic = zeros(1,size(maxlines,2));
breakpoints = [1 N+1];
breakpoints1 = breakpoints;
aic = AICcph(y,breakpoints,'poisson');
zvalues = zeros(1,size(maxlines,2));
[zvalue line] = max(abs(optimalcoef));
zvalues(line) = zvalue;
selectedlines = line;
optimalcoef(line) = 0;
remaininglines = find(abs(optimalcoef)>eps);
while ~isempty(remaininglines),
   bp = optimalposition(line);
   breakpoints = sort([breakpoints bp]);
   breakpoints1 = [breakpoints1 bp];
   aic = [aic AICcph(y,breakpoints,'poisson')];
   for rl = remaininglines,
      k = optimalposition(rl);
      if ismember(k,breakpoints),
         optimalcoef(rl) = 0;
      else
         l = optimalleft(rl);
         r = optimalright(rl);
         left = k-l+1;
         right = k+r;
         if k <  bp & right > bp,
            r = bp-k;
         end
         if k > bp & left < bp,
            l = k-bp+1;
         end
         [l,r,zu] = findnewoptimalunbalancedhaar(y,k,l,r,'poisson');
         optimalleft(rl) = l;
         optimalright(rl) = r;
         optimalcoef(rl) = zu;
      end
   end
   [zvalue line] = max(abs(optimalcoef));
   zvalues(line) = zvalue;
   optimalcoef(line) = 0;
   remaininglines = find(abs(optimalcoef)>eps);
   selectedlines = [selectedlines line];
end

A = length(aic);
daic = aic(2:A)-aic(1:A-1);
[minaic nCPaic] = min(aic);
nCPdaic = find(daic>-eps);
if isempty(nCPdaic),
   nCPdaic = length(aic);
else
   nCPdaic = nCPdaic(1);
                 % nPC is location of first minimum (i.e., location where daic
                 % is non-negative for the first time: daic is the difference
                 % between following and current aic value)
end
nCPdaic = nCPdaic-1;
nCPaic = nCPaic-1;
% nCPdaic and nCPaic are the number of change points to be added on top of the
% initial "change points" [1 N+1]

nCP = nCPaic;
cp = sort(breakpoints1(1:nCP+2));
for b = 1:nCP+1,
   bp0 = cp(b);
   bp1 = cp(b+1);
   ffaic(bp0:bp1-1) = mean(y(bp0:bp1-1));
end
cpaic = cp;

nCP = nCPdaic;
cp = sort(breakpoints1(1:nCP+2));
for b = 1:nCP+1,
   bp0 = cp(b);
   bp1 = cp(b+1);
   ffdaic(bp0:bp1-1) = mean(y(bp0:bp1-1));
end
cpdaic = cp;

stdevdaic = MAD1(daic)/0.6745;
outliers = find(daic<median(daic) - 3*stdevdaic);
cp = sort(breakpoints1([1 2 outliers+2]));
nCP = length(cp)-2;
for b = 1:nCP+1,
   bp0 = cp(b);
   bp1 = cp(b+1);
   ffoutliersdaic(bp0:bp1-1) = mean(y(bp0:bp1-1));
end
cpoutliersdaic = cp;
