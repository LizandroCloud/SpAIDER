if exist('gridnumber') ~= 1, gridnumber = 1; end
if exist('realrandom') ~= 1, realrandom = 0; end
if exist('samplesize') ~= 1, samplesize = 2067; end
N = samplesize;
x = makegrid(N,gridnumber,realrandom);

f = makesignal(x,'skyline');
y = makenoisy(f,1,0,10);
I = sparse(eye(N));
varerror = var2(y-f);
stdeverror = sqrt(varerror);
M = floor(log2(N));

% M = 5;

nt = 1;
n = 2;
typeofkernel = 'cos';
bandwidths = 6/N;


w = RWT_2GHtPU(x,y,M,'P','intpol','nt',nt);
Psitilde = RWT_2GHtPU(x,I,M,'P','intpol','nt',nt);

yy = IRT_2GHPU(x,w,M,'P','intpol','nt',nt);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','intpol','nt',nt);

maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Cubic interpolating prediction, no prefilter, 2-taps update

nt = 4;
n = 2;

w = RWT_2GHtPU(x,y,M,'P','intpol','nt',nt,'n',n);
Psitilde = RWT_2GHtPU(x,I,M,'P','intpol','nt',nt,'n',n);

yy = IRT_2GHPU(x,w,M,'P','intpol','nt',nt,'n',n);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','intpol','nt',nt,'n',n);

maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Cubic interpolating prediction, Unb.Haar prefilter, 2-taps update

w = RWT_2GHtPU(x,y,M,'P','intpol','nt',nt,'n',n,'Ht','UnbHaar');
Psitilde = RWT_2GHtPU(x,I,M,'P','intpol','nt',nt,'n',n,'Ht','UnbHaar');

yy = IRT_2GHPU(x,w,M,'P','intpol','nt',nt,'n',n);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','intpol','nt',nt,'n',n);

maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Next is Haar prediction (nt = 1) + Unb.Haar-prefilter; no update

w = RWT_2GHtPU(x,y,M,'P','intpol','nt',1,'Ht','UnbHaar');
Psitilde = RWT_2GHtPU(x,I,M,'P','intpol','nt',1,'Ht','UnbHaar');

yy = IRT_2GHPU(x,w,M,'P','intpol','nt',1);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','intpol','nt',1);

maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Kernel prediction, Unb.Haar-prefilter, no update

hJ = bandwidths;
w = RWT_2GHtPU(x,y,M,'P','kernel','type',typeofkernel,'h',hJ,'Ht','UnbHaar');
Psitilde = ...
     RWT_2GHtPU(x,I,M,'P','kernel','type',typeofkernel,'h',hJ,'Ht','UnbHaar');

yy = IRT_2GHPU(x,w,M,'P','kernel','type',typeofkernel,'h',hJ);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','kernel','type',typeofkernel,'h',hJ);
maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Local polynomial prediction, nt = 1; Unb.Haar-prefilter, no update
% Remark: current version of ThreshLab (3.2.1) has no final solution for local
% polynomial prediction when local neighborhood contains too few evens for
% non-singular local polynomial prediction; this results in NaN in the output.

ker = typeofkernel;
hJ = bandwidths;
nt = 1;

w = RWT_2GHtPU(x,y,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'Ht','UnbHaar');
Psitilde = ...
     RWT_2GHtPU(x,I,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'Ht','UnbHaar');

yy = IRT_2GHPU(x,w,M,'P','locpol','type',ker,'h',hJ,'nt',nt);
Iw = eye(length(w));
Psi = IRT_2GHPU(x,Iw,M,'P','locpol','type',ker,'h',hJ,'nt',nt);
maxerr =  max(abs(y-yy))
k = 26;
plot(x,Psi(1:N,k))

pause(2)

% Local polynomial prediction, nt = 1; biorthogonal pre-filter, no update
% next is very slow in current version (3.2.1)

w = RWT_2GHtPU(x,y,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'Ht','biortho');
Psitilde = ...
     RWT_2GHtPU(x,I,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'Ht','biortho');

Iw = eye(length(w));

y1 = IRT_2GHPU(x,w,M,'P','locpol','type',ker,'h',hJ,'nt',nt);
Psi1 = IRT_2GHPU(x,Iw,M,'P','locpol','type',ker,'h',hJ,'nt',nt);
maxerr1 =  max(abs(y-y1))

y2 = IRT_2GHPU(x,w,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'H','biortho');
Psi2 = IRT_2GHPU(x,Iw,M,'P','locpol','type',ker,'h',hJ,'nt',nt,'H','biortho');
maxerr2 =  max(abs(y-y2))


k = 56;
plot(x,Psi1(1:N,k))
plot(x,Psi2(1:N,k),'r')

pause(2)

