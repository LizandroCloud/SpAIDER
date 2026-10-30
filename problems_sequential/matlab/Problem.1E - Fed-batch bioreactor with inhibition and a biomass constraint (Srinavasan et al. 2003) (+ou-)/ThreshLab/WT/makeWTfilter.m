
function [h, ht] = makeWTfilter(type,n,nt);

%   Usage
%     [h ht] = makeWTfilter(biort-type,n,nt);
%     h = makeWTfilter(ort-type,n);
%   Inputs
%     ort-type   'Haar','Coiflet','Daubechies','Symmlet'
%     biort-type 'CDF', 'CDFless' (CDF with less dissimilar lengths)
%     n          number of (primal)  vanishing moments
%     nt         number of dual vanishing moments
%   Outputs
%     h          (primal) scaling filter
%     ht         dual scaling filter
%   Description
%     creates scaling function filter coefficients
%   Note
%     see help CDFspline for remark on precise definition of primal and dual
%     vanishing moments - your definition may just be the dual of the one
%     adopted here
%   See also
%     help FWT_BIO
%     help mirrorshifts
%     help CDFspline
%     help Haar
%     help Daubechies
%     help symmlet
%     help coiflet

if strcmp(type,'Haar') | strcmp(type,'haar'),
        [h ht] = Haar;
end

if strcmp(type,'Coiflet') | strcmp(type,'coiflet') | strcmp(type,'coif'),
        [h ht] = coiflet(n);
end

if strcmp(type,'Daubechies') | strcmp(type,'daubechies') | ...
   strcmp(type,'daub') | strcmp(type,'db') | strcmp(type,'Daub'),
        [h ht] = Daubechies(n);
end

if strcmp(type,'Symmlet') | strcmp(type,'symmlet') | strcmp(type,'symm') | ...
   strcmp(type,'Symlet') | strcmp(type,'symlet') | strcmp(type,'sym'),
   [h ht] = symmlet(n);
end

if strcmp(type,'CDF') | strcmp(type,'cdf') | strcmp(type,'CDFspline') | ...
   strcmp(type,'cdfspline'),
   [h ht] = CDFspline(n,nt);
end

if strcmp(type,'CDFLDL') | strcmp(type,'CDFless'),
   [h ht] = CDF_LDL(n,nt);
end

if nargin == 2,
   h = h/norm(h);
   ht = h;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
