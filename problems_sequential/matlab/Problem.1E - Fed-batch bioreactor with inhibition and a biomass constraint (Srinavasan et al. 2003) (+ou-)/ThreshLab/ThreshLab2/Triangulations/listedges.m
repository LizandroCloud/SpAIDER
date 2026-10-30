
function [E, etri, tE] = listedges(tri)

% listedges - ThreshLab2/Triangulations --- lists edges of triangulation
% Usage
%    E = listedges(tri);
%    [E etri] = listedges(tri);
%    [E, etri, tE] = listedges(tri);
% Inputs
%   tri   three-column matrix with indices (=vertices of triangulation)
% Outputs
%   E     list of edges, i.e., pairs of vertices. Vertex with lowest index 
%         always on the left hand side.
%   etri  list of indices into rows of tri, i.e. triangles in which the
%         corresponding edge in E occurs. etri has 2 columns, the element in 
%         the second column is zero if the corresponding edge lies on the 
%         boundary of the triangulation.
%   tE    list of indices into rows of E, i.e., a list of edges for every
%         triangle. tE has 3 columns with non-zero entries.
% 

ntri = size(tri,1);
edgelist = [tri(:,1:2)
         tri(:,2:3)
         tri(:,1:2:3)];
if isempty(edgelist),
   edgelist = zeros(0,4);
else
   edgelist = [edgelist repmat((1:ntri)',3,1) zeros(3*ntri,1)];
end
edgelist(:,1:2) = sort(edgelist(:,1:2),2); % sorts within each edge (vertex
                                           % with lowest index first)
edgelist = sortrows(edgelist); % sorts different edges
k = size(edgelist,1);
doubles = zeros(k,1);
doubles(2:k) = all(edgelist(1:k-1,1:2)'==edgelist(2:k,1:2)');
d = find(doubles);
r = find(1-doubles);
edgelist(d-1,4) = edgelist(d,3);
edgelist = edgelist(r,:);
doubles = (edgelist(:,1) == edgelist(:,2));
r = find(1-doubles);
edgelist = edgelist(r,:);
E = edgelist(:,1:2);
etri = edgelist(:,3:4);
tE = zeros(ntri,3);
tEi = ones(ntri,1);
if nargout > 2, for e = 1:length(E),
   t = edgelist(e,3);
   if t,
      tE(t,tEi(t)) = e;
      tEi(t) = tEi(t) + 1;
   end
   t = edgelist(e,4);
   if t,
      tE(t,tEi(t)) = e;
      tEi(t) = tEi(t) + 1;
   end
end, end
