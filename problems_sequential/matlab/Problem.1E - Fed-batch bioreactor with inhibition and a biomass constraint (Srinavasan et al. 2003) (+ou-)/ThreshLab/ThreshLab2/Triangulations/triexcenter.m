
function [C,r] = triexcenter(P1,P2,P3);

% triexcenter - ThreshLab2/Triangulation --- finds excenter of triangles
%   Usage
%      [C,r] = triexcenter(P1,P2,P3);
%      [C,r] = triexcenter(x3,y3);
%   Inputs
%      P1, P2, P3  matrices with two columns or two rows
%      x3,y3       matrices with three columns or three rows
%   Outputs
%      C   matrix of size equal to size(P1), size(P2), size(P3): center points
%      r   radiuses
%   Note
%      As the excenter of a triangle is also the center of circle in its
%      vertices, this routine can be seen as a vector version of circle3points
%   See also
%      help triincenter
%      help triorthocenter
%      help circle3points

transpose = 0;
if nargin == 3,
   sizeP1 = size(P1); sizeP2 = size(P2); sizeP3 = size(P3);
   if sizeP1~=sizeP2 | sizeP1~=sizeP3 | sizeP2~=sizeP3,
      error('sizes of input mismatch')
   end

   if sizeP1(2) == 2,
      P1 = P1'; P2 = P2'; P3 = P3'; transpose = 1;
   elseif sizeP1(1) ~= 2,
      error('threefold input must have two columns or two rows')
   end
   n = size(P1,2);

   x1 = P1(1,1:n);
   x2 = P2(1,1:n);
   x3 = P3(1,1:n);
   y1 = P1(2,1:n);
   y2 = P2(2,1:n);
   y3 = P3(2,1:n);
elseif nargin == 2,
   x3 = P1; y3 = P2;
   sizex3 = size(x3); sizey3 = size(y3);
   if sizex3 ~= sizey3,
      error('sizes of input mismatch')
   end

   if sizex3(2) == 3,
      x3 = x3'; y3 = y3'; transpose = 1;
   elseif sizex3(1) ~= 3,
      error('twofold input must have three columns or three rows')
   end
   n = size(x3,2);
   x1 = x3(1,1:n); x2 = x3(2,1:n); x3 = x3(3,1:n);
   y1 = y3(1,1:n); y2 = y3(2,1:n); y3 = y3(3,1:n);
end

xx = (x1+x2+x3)/3;
yy = (y1+y2+y3)/3;
x1 = x1-xx; x2 = x2-xx; x3 = x3-xx;
y1 = y1-yy; y2 = y2-yy; y3 = y3-yy;

a = x1.*y2+x2.*y3+x3.*y1-y1.*x2-y2.*x3-y3.*x1;
d = -((x1.^2+y1.^2).*y2+(x2.^2+y2.^2).*y3+(x3.^2+y3.^2).*y1) + ...
     ((x2.^2+y2.^2).*y1+(x3.^2+y3.^2).*y2+(x1.^2+y1.^2).*y3);
e =  ((x1.^2+y1.^2).*x2+(x2.^2+y2.^2).*x3+(x3.^2+y3.^2).*x1) - ...
     ((x2.^2+y2.^2).*x1+(x3.^2+y3.^2).*x2+(x1.^2+y1.^2).*x3);
f = -((x1.^2+y1.^2).*x2.*y3+(x2.^2+y2.^2).*x3.*y1+(x3.^2+y3.^2).*x1.*y2) + ...
     ((x1.^2+y1.^2).*y2.*x3+(x2.^2+y2.^2).*y3.*x1+(x3.^2+y3.^2).*y1.*x2);

o = find(abs(a)<eps);
a(o) = 1;
C = -[d./(2*a) ; e./(2*a)];
r = sqrt((d.^2+e.^2)/(4*a.^2) - f./a);
C(1:2,o) = NaN;
r(o) = Inf;
C = C + [xx; yy];
if transpose, C = C'; r = r'; end
