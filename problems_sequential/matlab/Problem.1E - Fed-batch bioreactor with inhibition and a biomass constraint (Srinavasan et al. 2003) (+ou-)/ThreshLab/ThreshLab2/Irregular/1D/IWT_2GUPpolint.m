
function y = IWT_2GUPpolint(x,w,nt,n,M);

% Inverse Wavelet Transform, 2nd Generation, update first

transposex = 0;
if size(x,2) > size(x,1),
   x = x';
end
transposew = 0;
if size(w,1) ~= size(x,1),
   w = w';
end

moments = cell(M);
grids = cell(M);
if n > 0,
   Mfine = initmom(x,(0:n-1));
else
   Mfine = zeros(length(x),0);
end
xfine = x;

for j = 1:M
   moments{j} = Mfine;
   grids{j} = xfine;
   [even, odd] = splitevenodd(xfine);
   xeven = xfine(even);
   xodd = xfine(odd);
   Ne = length(even);
   Ie = sparse(eye(Ne));
   Meven = Mfine(even,1:n); Modd = Mfine(odd,1:n);
   e1 = ones(size(xeven)); o1 = ones(size(xodd));
   rescale = e1 + updateprimalmoments(xodd,o1,xeven,Modd,Meven);
   D = sparse(diag(rescale));
   Mcoarse = D*Meven;
   Mfine = Mcoarse;
   xfine = xeven;
end

y = w;
yeven = y(1:Ne,:);
for j = M:-1:1,
   xfine = grids{j};
   [even, odd] = splitevenodd(xfine);
   xeven = xfine(even);
   xodd = xfine(odd);
   Mfine = moments{j};
   Nj = length(xfine);
   Ne = length(xeven);
   No = Nj-Ne;

   yodd = w(Ne+1:Nj,:);
   Ie = sparse(eye(Ne));
   P = sparse(polynomialprediction(xeven,Ie,xodd,nt));
   % PREDICTION %
   yodd = yodd + P*yeven;
   % UPDATE
   Modd = Mfine(odd,:);
   Meven = Mfine(even,:);
   if n > 0,
      e1 = ones(size(yeven)); o1 = ones(size(yodd));
      rescale = e1 + updateprimalmoments(xodd,o1,xeven,Modd,Meven);
      yeven = yeven.*rescale;
      yeven = yeven - updateprimalmoments(xodd,yodd,xeven,Modd,Meven);
   end
   % MERGE %
   y(even,:) = yeven;
   y(odd,:) = yodd;
   yeven = y(1:Nj,:);
end

if transposew, y = y'; end
