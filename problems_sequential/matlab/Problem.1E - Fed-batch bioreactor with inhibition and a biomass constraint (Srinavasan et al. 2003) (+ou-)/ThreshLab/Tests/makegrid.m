function x = makegrid(N,number,realrandom)

% usage: x = makegrid(N,number,realrandom)
% the argument realrandom is optional. Its default value is 0.

if nargin < 3 | ~realrandom,
   rand('seed',931316785);
end

if nargin < 2, number = 0; end

switch number
case 0
   x = (1:N)/N;
case 1
   x = sort(rand(1,N));  
case 2
   N1 = floor(N/20);
   N2 = N1 + floor(N/200);
   x(1:N1) = rand(size(1:N1));
   x(1:N1) = x(1:N1)/5;
   x(N1+1:N2) = rand(size(N1+1:N2));
   x(N1+1:N2) = 0.2 + x(N1+1:N2)*0.2;
   x(N2+1:N) = rand(size(N2+1:N));
   x(N2+1:N) = 0.4 + x(N2+1:N)*0.6;
   x = sort(x);
case 3
   N1 = floor(N/20);
   N2 = N1 + floor(N/100);
   x(1:N1) = rand(size(1:N1));
   x(1:N1) = x(1:N1)*0.2;
   x(N1+1:N2) = rand(size(N1+1:N2));
   x(N1+1:N2) = 0.2 + x(N1+1:N2)*0.2;
   x(N2+1:N) = rand(size(N2+1:N));
   x(N2+1:N) = 0.4 + x(N2+1:N)*0.6;
   x = sort(x);
case 31
% sparsely sampled interval in the middle
   N1 = floor(47*N/100);
   N2 = N1 + floor(N/20);
   N3 = N2 + floor(N/100);
   x(1:N1) = rand(size(1:N1));
   x(1:N1) = x(1:N1)*0.3;
   x(N1+1:N2) = rand(size(N1+1:N2));
   x(N1+1:N2) = 0.3 + x(N1+1:N2)*0.2;
   x(N2+1:N3) = rand(size(N2+1:N3));
   x(N2+1:N3) = 0.5 + x(N2+1:N3)*0.2;
   x(N3+1:N) = rand(size(N3+1:N));
   x(N3+1:N) = 0.7 + x(N3+1:N)*0.3;
   x = sort(x);
case 32
   % normally distributed samples
   x = sort(randn(1,N));
   x = (x-min(x))/(max(x)-min(x));
case 4
   N3 = floor(N/20);
   N2 = floor(N/100);
   N1 = floor(N/2) - N2 - N3;
   N2 = N1+N2;
   N3 = N2+N3;
   x(1:N1) = rand(size(1:N1));
   x(1:N1) = x(1:N1)/4;
   x(N1+1:N2) = rand(size(N1+1:N2));
   x(N1+1:N2) = 0.25 + x(N1+1:N2)*0.25;
   x(N2+1:N3) = rand(size(N2+1:N3));
   x(N2+1:N3) = 0.5 + x(N2+1:N3)*0.25;
   x(N3+1:N) = rand(size(N3+1:N));
   x(N3+1:N) = 0.75 + x(N3+1:N)*0.25;
   x = sort(x);
case 5 % sort of inverse devil's staircase
   k = N/(N-1);
   x = (0:1);
   nn = length(x);
   while (nn < N)
      xeven = x;
      xodd = xeven(1:nn-1);
      xodd = xodd + rand(size(xodd)).*(xeven(2:nn)-xeven(1:nn-1))/k;
      clear x
      nn = 2*nn-1;
      x(1:2:nn) = xeven;
      x(2:2:nn) = xodd;
      xeven = x;
      xodd = xeven(2:nn);
      xodd = xodd - rand(size(xodd)).*(xeven(2:nn)-xeven(1:nn-1))/k;
      clear x
      nn = 2*nn-1;
      x(1:2:nn) = xeven;
      x(2:2:nn) = xodd;
   end
   x = x(1:N)/x(N);
case 50 % non-random version of case 5
   k = (3/2)*N/(N-1);
   % k = 3*N/(N-1);
   x = (0:1);
   nn = length(x);
   while (nn < N)
      xeven = x;
      xodd = xeven(1:nn-1) + (xeven(2:nn)-xeven(1:nn-1))/k;
      clear x
      nn = 2*nn-1;
      x(1:2:nn) = xeven;
      x(2:2:nn) = xodd;
      xeven = x;
      xodd = xeven(2:nn) - (xeven(2:nn)-xeven(1:nn-1))/k;
      clear x
      nn = 2*nn-1;
      x(1:2:nn) = xeven;
      x(2:2:nn) = xodd;
   end
   x = x(1:N)/x(N);
case 8
   for (i = 1:9) x(i) = (i-1)/8; end
   M = 9;
   while (N - M > 3)
      for (i = 1:4) x(M+i) = (x(i) + x(i+1))/2; end
      M = M + 4;
      x = sort(x);
   end  
case 9
   Gap = 1000;
   k = div(N,4)*2+1;
   x = (1:N);
   x(k) = x(k) + Gap;
   x(k+1:N) = x(k+1:N) + 2*Gap;
   x = x/x(N);
case 610 % extreme case of grid with given homogeneity constant gamma>1; 
         % See Daubechies e.a., Regularity of irregular subdivision (Constr.
	 % Approx.)
   gamma = 5;
   x = makextrmhomogeneousgrid(N,gamma);
case 620 % non-homogeneous, but dyadically balanced grid with given balance
         % beta
   beta = 1/6;
   % x = (0:1);
   x = (0:2)/2;
   nn = length(x);
   while nn < N,
      xx(1:2:2*nn-1) = x;
      xx(2:2:2*nn-1) = beta*x(1:nn-1) + (1-beta)*x(2:nn);
      nn = length(xx);
      x = xx;
   end
   x = x(1:N)/x(N);
case 666 % devil's staircase, cantor function
   x = makecantorgrid(N);
case 6660 % a generalisation of case 666 (666 is special case with Gap = 1/3)
   Gap = 1/10;
   x = makecantorgrid(N,Gap);
case 6661 % random version of case 666
   Gap = 1/10;
   x = makecantorgrid(N,Gap,1,realrandom);
otherwise
  x = (1:N)/N;
end
