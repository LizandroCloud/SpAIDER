
function B = initboundary(x,y,anglethreshold)

% initboundary - ThreshLab2/Triangulations --- initialise triangulation
%   Usage
%     B = initboundary(x,y,anglethreshold);
%   Inputs
%     x     x-coordinates of 2d-scattered data
%     y     y-coordinates of 2d-scattered data
%     length(x) = length(y)
%     anglethreshold optional argument; default = pi/2 (see below)
%   Outputs
%     B     ordered list of indices into (x,y), the vertices on the
%           boundary of the scattered data
%   Description
%     constructs a boundary of data cloud (x,y) such that no triangles can be
%     constructed connecting two adjacent boundary points with an inside point,
%     where the inside angle is larger than anglethreshold.
%     This routine does not use Delaunay or any other triangulation
%   See also
%     help inittri
%     help removenarrowtrinearboundary
%     help boundary
%     help convexhull2

x = row(x); y = row(y);
N = length(x);
if length(y) ~= N, error('dimension mismatch x and y'); end

if nargin < 3, anglethreshold = pi/2; end
cosanglethreshold = cos(anglethreshold);

B = convexhull2(x,y);
nB = length(B);
if B(1) ~= B(nB), B = [B B(1)]; end
anglesB = anglesboundary(B,x,y);

if length(anglesB) == length(B)-1, anglesB = [anglesB anglesB(1)]; end

nB = length(B);
C = ones(1,nB-1);
while any(C),
   c = min(find(C));
   p1 = B(c); p2 = B(c+1); pmin = min(p1,p2); pmax = max(p1,p2);
   % p0 = [(1:pmin-1) NaN (pmin+1:pmax-1) NaN (pmax+1:N)];
   p0 = (1:N);
   x1 = x(p1)-x(p0); y1 = y(p1)-y(p0);
   x2 = x(p2)-x(p0); y2 = y(p2)-y(p0);
   p1p2 = ((x1.*x2)+(y1.*y2))./(sqrt(x1.^2+y1.^2).*sqrt(x2.^2+y2.^2));
   p1p2(B) = NaN; % existing points on boundary excluded from re-use
   [cosmaxangle p0] = min(p1p2);
   if cosmaxangle < cosanglethreshold,
      % anglesB = anglesboundary(B,x,y);
      anglesp1p0p2 = anglesboundary([p1 p0 p2],x,y);
      if anglesB(c) > anglesp1p0p2(1) & anglesB(c+1) > anglesp1p0p2(3),
         % the removed triangle must fit inside the current boundary, otherwise
         % the new boundary will have self-intersections
         C = [C(1:c) 1 C(c+1:nB-1)];
         B = [B(1:c) p0 B(c+1:nB)];
         anglesB = [anglesB(1:c) 2*pi anglesB(c+1:nB)];
         anglesB(c:c+2) = anglesB(c:c+2) - anglesp1p0p2;
         nB = nB+1;
      else
         C(c) = 0;
      end
   else
      C(c) = 0;
   end
end
