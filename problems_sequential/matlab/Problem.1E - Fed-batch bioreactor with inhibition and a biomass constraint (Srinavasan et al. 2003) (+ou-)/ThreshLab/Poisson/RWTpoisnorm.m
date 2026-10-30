
function [pn, firstw, lastw] = RWTpoisnorm(x,M,ht,gt,sht,sgt)

% usage: pn = RWTpoisnorm(x,M,ht,gt,sht,sgt);
% or: [pn firstw lastw] = RWTpoisnorm(x,M,ht,gt,sht,sgt);
% computes appropriate normalisation for wavelet coefficients of Poisson noise
% vector x. The parameters L,ht,gt,sht,sgt correspond to those in the function
% call wcoef = RWT_BIO(x,L,ht,gt,sht,sgt);
% firstw is index into first sample in x participating in the normalisation
% factor; similar for lastw

N = length(x);
nht = length(ht);
ngt = length(gt);
firsts = (1:N);
lasts = firsts;
firstw = zeros(M+1,N);
lastw = zeros(M+1,N);
J = ceil(log2(N));
for j = 1:M,
   firstw(j,:) = firsts(mod((1:N)+sgt-1,N)+1);
   lastw(j,:)  = lasts(mod((1:N)+ngt+sgt-2,N)+1);
   firsts = firsts(mod((1:N)+sht-1,N)+1);
   lasts  = lasts(mod((1:N)+nht+sht-2,N)+1);
   sht = sht*2;
   sgt = sgt*2;
   nht = nht*2 - 1;
   ngt = ngt*2 - 1;
end
firstw(M+1,:) = firsts;
lastw(M+1,:) = lasts;
for j = 1:M+1,
   support = length(WTbasisfunction(j,ht,gt,sht,sgt,0));
   if support >= N, % If the basis functions have support that covers the whole
                    % interval, the calculation with first (begin point) and
                    % last (end point) fails, since the algorithm then only
                    % considers the observations between those points, whereas
                    % the right value of the normalisation factor involves all
                    % the observations of course.
       pn(j,:) = sum(x); 
   else
      for k = 1:N,
         if firstw(j,k) > lastw(j,k),
            pn(j,k) = sum(x(firstw(j,k):N)) + sum(x(1:lastw(j,k)));
         else
            pn(j,k) = sum(x(firstw(j,k):lastw(j,k)));
         end
      end
   end
end
pn = sqrt(abs(pn));

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
