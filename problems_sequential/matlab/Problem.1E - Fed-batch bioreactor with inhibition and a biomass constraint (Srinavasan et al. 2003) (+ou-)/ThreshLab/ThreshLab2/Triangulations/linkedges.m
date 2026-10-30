
function [points e] = linkedges(E);

% linkedges - ThreshLab2/Triangulations --- link edges in correct order
%   Usage
%     [points e] = linkedges(E);
%   Inputs
%     E  two-column vector with indices = vertices of edges to be linked
%   Outputs
%     points  ordered list of linked points (= polyline)
%     e       index into rows of E, indicating order of selection of edges in E
%   Description
%     -  E is supposed to be "linkable", otherwise, gaps are detected and 
%     reported as error message
%     -  linkedges can deal with self-intersections in polylines
%     -  closed polylines have points(1) == points(last)
%     -  Unless the edges do not constitute a closed loop AND the first edge in
%     E does not contain an end point of the chain, the list in e ALAWYS starts
%     with a 1. (and the list of points starts accordingly, of course)
%   See also
%     

n = size(E,1);
e = (1:n);
points = zeros(1,n+1);
if ~isempty(E),
   points(1) = E(1,1);
   points(2) = E(1,2);
   E = E(2:n,:);
end
for k = 2:n
   rr = find(sum(E==points(k),2));
   if isempty(rr)
      e(1:k-1) = e(k-1:-1:1);
      points(1:k) = points(k:-1:1);
      rr = find(sum(E==points(k),2));
   end
   last = points(k);
   nshifts = 0;
   while isempty(rr), % if rr is till empty after flipping, this means that we
                      % have ended up in a closed loop. We have to look for the
                      % exit
      points(1:k-1) = shift(points(1:k-1));
      e(1:k-1) = shift(e(1:k-1));
      rr = find(sum(E==points(1),2));
      last = points(1);
      nshifts = nshifts+1;
      if nshifts > k, % if we have tried all shifts, this means that the edges
                      % do not connect all vertices. There is a gap.
                      % create an artificial edge bridging the gap.
       error('linkedges has detected a gap, edges do not connect all vertices')
      end
   end
   e(k) =  rr;
   % This index in e is not yet (necessarily) correct: it currently carries a 
   % relative number, dependent on the state of E after elimination of 
   % previously selected edges, we now reconstruct the absolute order.
   ee = e(1:k-1);
   l = find(ee<=e(k));
   while ~isempty(l),
      e(k) = e(k)+1;
      ee = setdiff(ee,ee(l(1)));
      l = find(ee<=e(k));
   end
   points(k) = last;
   if (E(rr(1),1) == points(k))
      points(k+1) = E(rr(1),2);
   else
      points(k+1) = E(rr(1),1);
   end
   E = [E(1:rr(1)-1,:)' E(rr(1)+1:n-k+1,:)']';
end
if points(n+1) == points(1)
   points = points(1:n);
end
