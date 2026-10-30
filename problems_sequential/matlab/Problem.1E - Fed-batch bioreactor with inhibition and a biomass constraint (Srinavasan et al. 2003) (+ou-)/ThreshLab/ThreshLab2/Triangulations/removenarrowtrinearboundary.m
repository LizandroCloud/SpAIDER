
function [tri, B] = removenarrowtrinearboundary(tri,x,y,B,anglethreshold)

% removenarrowtrinearboundary -- ThreshLab2/Triangulations
%   Usage
%      [tri, B] = removenarrowtrinearboundary(tri,x,y,B,anglethreshold);
%      tri = removenarrowtrinearboundary(tri,x,y,B);
%      tri = removenarrowtrinearboundary(tri,x,y);
%      tri = removenarrowtrinearboundary(tri,x,y,anglethreshold);
%   Inputs
%      tri   triangulation, i.e., triple indices into (x,y)
%      x,y   coordinates of vertices
%      B     boundary or section of boundary, i.e., vector of vertices 
%            (indices into x,y) along which removal takes place
%            default is output of triboundary(tri)
%      anglethreshold  maximum value of angle in triangle on boundary 
%                      (default pi/2)
%   Outputs
%      points  triangulation after removal
%      B       boundary after removal
%   See also
%      help boundary
%      help inittri
%   Notes
%      if B represents a closed loop, B(last) MUST be B(1) !!

x = row(x);
y = row(y);

if nargin == 4,
   if length(B) < 2,
      anglethreshold = B;
      B = triboundary(tri);
      B = [B B(1)];
   else
      anglethreshold = pi/2;
   end
end
if nargin == 3,
   B = triboundary(tri);
   B = [B B(1)];
   anglethreshold = pi/2;
end

if length(B)<2,
   error('Boundary must be a vector: help removenarrowtrinearboundary')
end

ntri = size(tri,1);
nB = length(B);
nE = nB - 1;
E = [column(B(1:nE)) column(B(2:nE+1))];
indexE1inB = (1:nE);
while nE,
   facet = E(1,:);
   t = ismember(tri,facet);
   t = sum(t');
   t = (t==2);
   t = find(t);
   if length(t) > 1,
      error('found more than one triangle on boundary facet')
   elseif length(t) == 1,
      narrow = (anglesboundary(tri(t,:),x,y) > anglethreshold);
      thirdpoint = setdiff(tri(t,:),facet);
      inside = ~ismember(thirdpoint,B);
      thirdpos = ~ismember(tri(t,:),facet);
      if any(and(thirdpos,narrow)) & inside,
         tt = setdiff((1:ntri),t);
         tri = tri(tt,:);
         ntri = ntri - 1;
         Enew = [facet(1) thirdpoint
                 thirdpoint facet(2)];
         E = [E' Enew']';
         nE = nE+2;

         b = indexE1inB(1);
         B = [B(1:b) thirdpoint B(b+1:nB)];
         indexE1inB(indexE1inB>b) = indexE1inB(indexE1inB>b) + 1;
         indexE1inB = [indexE1inB b b+1];
         nB = nB + 1;
      end
   end
   E = E(2:nE,:);
   indexE1inB = indexE1inB(2:nE);
   nE = nE-1;
end
