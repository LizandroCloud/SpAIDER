function tri = bndrycstrbndrydelaunay1(points,x,y);

% bndrycstrbndrydelaunay1 - ThreshLab2/Triangulations ---
%                       boundary constrained boundary delaunay triangulation
%   Usage
%      tri = bndrycstrbndrydelaunay1(points,x,y);
%   Inputs
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%      points     ordered list of indices into x and y corresponding to vertices
%   Outputs
%      tri        triangulation
%   Description
%      triangulation within boundary described by points. All vertices of
%      points are assumed to be on the boundary (if this is not the case, use
%      bndrycstrdelaunay).
%
%      The basic idea is that we know that there must be a triple of
%      successive vertices in the list of points that forms a triangle of the
%      eventual delaunay triangulation. This procedure looks for such a set and
%      then applies the same trick recursively to the set of points after
%      removing the central point of the just found triangle from the list of
%      points.
%   See also
%      help bndrycstrbndrydelaunay
%      help bndrycstrdelaunay1


tri = zeros(0,3);
n = length(points);
while n > 3,
   pp = [points(n) points points(1)];
   alfa = anglesboundary(points,x,y);
   p = 0;
   trianglefound = 0;
   while ~trianglefound,
      p = p+1;
      if alfa(p) < pi,
         possibletriangle = pp(p:p+2);
         P1 = [x(pp(p))     y(pp(p))];
         P2 = [x(pp(p+1)) y(pp(p+1))];
         P3 = [x(pp(p+2)) y(pp(p+2))];
         otherpoints = setdiff(points,possibletriangle);
         [C,r] = circle3points(P1,P2,P3);
         % [C,r] = triexcenter(x(pp(p:p+2)),y(pp(p:p+2)));
         x0 = C(1);
         y0 = C(2);
         otherondisk = pointondisk(x(otherpoints),y(otherpoints),x0,y0,r);
         for q = 1:length(otherpoints),
            % any other point on the disk circumscribing the possible triangle
            % is eliminated from the list of obstacles if the edge connecting 
            % that point to the opposite point in the candidate triangle is not
            % within the region to be triangulated
            if otherondisk(q),
               poq = [points(p) otherpoints(q)];
               otherondisk(q) = edginside(poq,points,x,y);
            end
         end
         if ~any(otherondisk),
            trianglefound = 1;,
            tri = [tri' possibletriangle']';
            points = [points(1:p-1) points(p+1:n)];
            n = n-1;
         end
      end
   end
end
tri = [tri' points']';
