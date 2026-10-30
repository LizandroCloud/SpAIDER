function grey = grey(n);

% usage: grey = grey(n);
% creates matlab colormap for grey value images with n levels;
% default for n = 256.

if nargin<1,
   n = 256;
end

grey = [(0:n-1)' (0:n-1)' (0:n-1)']/n;
