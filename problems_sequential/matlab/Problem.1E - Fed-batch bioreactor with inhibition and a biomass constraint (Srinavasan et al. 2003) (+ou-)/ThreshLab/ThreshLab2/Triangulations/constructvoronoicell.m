
function [voronoicell neighbours] = constructvoronoicell(v0,x,y)

% constructvoronoicell - ThreshLab2/Triangulations --- Voronoi cell for single
%                                                    vertex
%   Usage
%     [v n] = constructvoronoicell(v0,x,y)
%   Inputs
%     v0    vertex; index into vectors x and y
%     x     x-coordinates of 2d-scattered data
%     y     y-coordinates of 2d-scattered data
%     length(x) = length(y)
%   Outputs
%     vorornoicell  two column matrix; each row containing the coordinates of
%                   one corner of the voronoi cell of v0
%     neighbours    ordered list of indices into x and y containing the
%                   adjacent vertices (i.e., vertices whose voronoi cell
%                   borders the voronoi cell of v0)
%   Description
%     vertices in x,y are ordered according to distance to v0
%     new neighbours are investigated as long as there are vertices at distance
%     from v0 smaller than twice the furthest point in the current voronoi cell
%   See also

x = row(x); y = row(y);
npoints = length(x);
if length(y) ~= npoints, 
   warning('vectors x and y must have equal length. Ignoring abundant points')
   npoints = min(length(x),length(y));
   x = x(1:npoints);
   y = y(1:npoints);
end

B = convexhull2(x,y);
xB = x(B); yB = y(B); XB = [xB;yB];
YB = shiftsegmentspolygone(XB,1.1);
xB = YB(1,:); yB = YB(2,:);
x = [x xB];
y = [y yB];
% All points in the original (x,y) are now strictly inside the convex hull of
% the new (x,y), and thus their voronoi cells cannot be unbounded. Further
% calculations assume boundedness of the cells.
npoints = length(x);

x0 = x(v0); y0 = y(v0);
d2 = (x-x0).^2 + (y-y0).^2;
[d2 closestpoints] = sort(d2);
furthestcorner = Inf;
k = 1;
bisectors = zeros(0,3); 
% rows of bisectors contain coefficients a b c in equation ax+by+c = 0
corners = zeros(0,5);
% rows of corners contain x-coordinate, y-coordinate, two indices into
% bisectors indicating the lines crossing in that corner point  and a binary
% variable which is one if that corner point is active, i.e., still on the
% corner of the current voronoi cell
while furthestcorner > d2(k+1)/4,
% stop criterion: distance to any other remaining vertex is larger than twice
% the distance to most remote point inside current cell. This may last a while 
% for points on boundary with large voronoi cells.
   n = closestpoints(k+1);
   xn = x(n); yn = y(n);
   d  = [xn-x0 yn-y0]; % -d(1)/d(2) is slope of perpendicular bisector
   xm = (xn+x0)/2; ym = (yn+y0)/2; % midpoint
   abc = [d(1) d(2) -xm*d(1)-ym*d(2)];
           % abc are coefficients of line ax+by+c = 0
   abc = abc * sign(abc*[x0;y0;1]);
   if size(corners,1) > 0,
      S = [corners(:,1:2) ones(size(corners,1),1)]*abc';
      corners(:,5) = corners(:,5) & (S>-eps);
   end
   for b = 1:k-1,
      A = [abc(1:2);bisectors(b,1:2)]; 
      if cond(A) < 10^10,
         R = -[abc(3);bisectors(b,3)]; C = A\R;
         S = bisectors * [C; 1];
         corners = [corners; C(1) C(2) b k all(S>-eps)];
      end
   end
   k = k+1;
   bisectors = [bisectors; abc];
   ac = find(corners(:,5)); nac = length(ac);
   if nac > 2,
      % activecorners = corners(ac,1:4);
      xc = corners(ac,1); yc = corners(ac,2);
      dc2 = (xc-x0).^2 + (yc-y0).^2;
      furthestcorner = max(dc2);
      neighbours = corners(ac,3:4);
      neighbours = closestpoints(neighbours+1);
      neighboursrow = sort(reshape(neighbours,1,2*nac));
      if any(neighboursrow(1:2:2*nac)~=neighboursrow(2:2:2*nac)),
         furthestcorner = Inf;
      else
         [neighbours i] = linkedges(neighbours);
      end
   end
   if k > npoints-2, furthestcorner = 0; end
end
voronoicell = corners(ac,1:2);
voronoicell = voronoicell(i,1:2);
