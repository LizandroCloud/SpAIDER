
function [k,l,r,z] = findoptimalunbalancedhaar(y,k0,kmin,kmax,scale,pdf);

% findoptimalunbalancedhaar - ThreshLab2/ContinuousUnbalancedHaar
%   Usage
%     [k,l,r,z] = findoptimalunbalancedhaar(y,k0,kmin,kmax,scale,pdf);
%   Inputs
%     y       observed (noisy) values
%     k0      location around which unbalanced analysis will be performed
%     kmin    left most location to be considered (e.g. given another change
%     kmax    right most location to be considered                     point)
%     scale   initial scale of balanced analysis (starting value for output
%             variables k and l)
%     pdf     string: 'normal' or 'poisson' for Haar or Haar-Fisz analysis
%             respectively
%   Outputs
%     k       center location of unbalanced Haar basis
%     l       left scale of unbalanced Haar basis
%     r       right scale of unbalanced Haar basis
%     z       signficance/Haar(or Haar-Fisz) wavelet coefficient
%   Description
%     finds optimal values for k,l and r for an unbalanced Haar or
%     Haar-Fisz transform:
%     (mean(y(k+1:k+r))-mean(y(k-l+1:k)))/sqrt(nn)*sqrt(l*r);
%     where nn = sum(y(k-l+1:k+r)); in the case of Haar-Fisz
%     and   nn = r+l                in the case of Haar
%     takes a maximum in the neighbourhood of k0. k0 may fluctuate (into output
%     k)
%   See also
%     help findallunbalancedhaar (calling routine)

N = length(y);
k = k0;
l = min(scale,k);
r = min(scale,N-k);
maxzu = 0;
maxl = l;
maxr = r;
lmin = 1;
lmax = k;
rmin = 1;
rmax = N-k;
lstep = ceil((lmax-lmin)/512);
rstep = ceil((rmax-rmin)/512);

while lstep > 0 | rstep > 0,
% if one of the steps equals zero, matlab will not go through the iteration (so
% there is no danger for inifinte loops from this point of view)
%lstep
%rstep

   for l = lmin:lstep:lmax,
      if strcmp(pdf,'poisson') | strcmp(pdf,'Poisson') | strcmp(pdf,'pois'),
         nn = sum(y(k-l+1:k+r));
      else
         nn = r+l;
      end
      if nn > 0 & r*l >0,
         zu = (mean(y(k+1:k+r))-mean(y(k-l+1:k)))/sqrt(nn)*sqrt(l*r);
         if abs(zu) > abs(maxzu),
            maxl = l;
            maxzu = zu;
         end
      end
   end
   l = maxl;
   for r = rmin:rstep:rmax,
      if strcmp(pdf,'poisson') | strcmp(pdf,'Poisson') | strcmp(pdf,'pois'),
         nn = sum(y(k-l+1:k+r));
      else
         nn = r+l;
      end
      if nn > 0 & r*l >0,
         zu = (mean(y(k+1:k+r))-mean(y(k-l+1:k)))/sqrt(nn)*sqrt(l*r);
         if abs(zu) > abs(maxzu),
            maxr = r;
            maxzu = zu;
         end
      end
   end
   r = maxr;
   for l = lmin:lstep:lmax,
      if strcmp(pdf,'poisson') | strcmp(pdf,'Poisson') | strcmp(pdf,'pois'),
         nn = sum(y(k-l+1:k+r));
      else
         nn = r+l;
      end
      if nn > 0 & r*l >0,
         zu = (mean(y(k+1:k+r))-mean(y(k-l+1:k)))/sqrt(nn)*sqrt(l*r);
         if abs(zu) > abs(maxzu),
            maxl = l;
            maxzu = zu;
         end
      end
   end
   l = maxl;
   lmin = max(l-lstep,lmin);
   lmax = min(l+lstep,lmax);
   rmin = max(r-rstep,rmin);
   rmax = min(r+rstep,rmax);
   lstep = max(0,min(lstep-1,ceil((lmax-lmin)/512)));
   rstep = max(0,min(rstep-1,ceil((rmax-rmin)/512)));

end

% now adjust/optimise central position within given left and right boundary
% this may hide second change points within close distance

left = k-l+1;
right = k+r;
kmin = max(kmin,left);
kmax = min(kmax,right);
kstep = ceil((kmax-kmin)/512);
while kstep > 0,

   maxk = k;
   left = k-l+1;
   right = k+r;
   if strcmp(pdf,'poisson') | strcmp(pdf,'Poisson') | strcmp(pdf,'pois'),
      nn = sum(y(left:right));
   else
      nn = right - left + 1;
   end
   for k = kmin:kstep:kmax,
      l = k-left+1;
      r = right-k;
      if l > 0 & r > 0,
         zu = (mean(y(k+1:right))-mean(y(left:k)))/sqrt(nn)*sqrt(l*r);
      else
         zu = 0;
      end
      if abs(zu) > abs(maxzu),
         maxk = k;
         maxl = l;
         maxr = r;
         maxzu = zu;
      end
   end
   k = maxk;
   l = maxl;
   r = maxr;

   kmin = max(k-kstep,kmin);
   kmax = min(k+kstep,kmax);
   kstep = min(kstep-1,ceil((kmax-kmin)/512));

end
z = maxzu;
