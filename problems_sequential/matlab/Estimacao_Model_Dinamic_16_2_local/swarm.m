function [xo,Ot,nS]=swarm(S,x0,ip,isb,Lb,Ub,R,nger,npas,problem,c1,c2)
%   Hard constrained global optimization using SWARM.
%
%   [xo,Ot,nS]=swarm(S,x0,ip,isb,Lb,Ub,R,nger,npas,problem,c1,c2)
%
%   S: objective function
%   x0: inital point
%   ip: (0): no plot (default), (>0) plot figure ip with pause, (<0) plot figure ip
%   isb: isb(i)=1 is i-th variable is bounded, otherwise isb(i)=0 (default = 0)
%   Lb, Ub: lower and upper bound vectors (default = x0*(1+/-2) for plot if isb(i)=0)
%   R: axis vector of the hyperellipse centered in x0 (default = max(0.1*abs(x0+~x0),1))
%   nger: number of generations or flights (default = 50)
%   npas: number of birds (default = 10)
%   problem: (-1): minimum (default), (1): maximum
%   c1: local search parameter (default = 1)
%   c2: global search parameter (default = 1)
%   xo: optimal point
%   Ot: optimal values of S
%   nS: number of objective function evaluations

%   $Revision: 1.0 $  $Date: 2002/09/15 23:10:05 $
%   Author: Esly Ferreira da Costa Junior (versao FORTRAN, 10/2001)
%   MATLAB version: Argimiro R. Secchi (arge@enq.ufrgs.br)

 if nargin < 2,
   error('swarm requires two input arguments');
 end

 x0=x0(:);
 n=size(x0,1);
 
 if nargin < 3 | isempty(ip),
   ip=0;
 end
 if nargin < 4 | isempty(isb),
   isb=zeros(n,1);
 end
 if nargin < 5 | isempty(Lb),
   Lb=-x0-~x0;
 end
 if nargin < 6 | isempty(Ub),
   Ub=2*x0+~x0;
 end
 if nargin < 7 | isempty(R),
   R=max(0.1*abs(x0+~x0),1);
 end
 if nargin < 8 | isempty(nger),
   nger=50;
 end
 if nargin < 9 | isempty(npas),
   npas=10;
 end
 if nargin < 10 | isempty(problem),
   problem=-1;
 end
 if nargin < 11 | isempty(c1),
   c1=1;
 end
 if nargin < 12 | isempty(c2),
   c2=1;
 end

% Local variables:
%   ipg  = index of best position
%   tw   = inertial factor reduction rate
%   p    = matrix of best position of each bird
%   v    = bird velocity matrix
%   x    = bird location matrix
%   w    = inertial factor

                    % initialization 
 R=abs(R(:));
 Lb=Lb(:);
 Ub=Ub(:);
 w=0.9;
 tw=(w-4e-3)/nger;
 Ot=feval(S,x0)*problem;
 y(1)=Ot;
 p=zeros(n,npas);
 v=zeros(n,npas);
 x=zeros(n,npas);
   
 if ip & n == 2,
   figure(abs(ip));
   [X1,X2]=meshgrid(Lb(1):(Ub(1)-Lb(1))/20:Ub(1),Lb(2):(Ub(2)-Lb(2))/20:Ub(2));
   [n1,n2]=size(X1);
   f=zeros(n1,n2);
   for i=1:n1,
    for j=1:n2,
      f(i,j)=feval(S,[X1(i,j);X2(i,j)]);
    end
   end
   mxf=max(max(f));
   mnf=min(min(f));
   df=mnf+(mxf-mnf)*(2.^(([0:10]/10).^2)-1);
   [vl,h]=contour(X1,X2,f,df); hold on;
%   clabel(vl,h);
   h1=plot(x0(1),x0(2),'ro');
   legend(h1,'start point');

   if ip > 0,
     disp('Pause: hit any key to continue...'); pause;
   end
 end
 
 del=Ub-Lb;
 
 for j=1:npas,
    x(:,j)=Lb+rand(n,1).*del;
    v(:,j)=0.1*rand(n,1).*del;
    p(:,j)=x(:,j);
 end	

 x(:,1)=x0;
 p(:,1)=x0;
 nS=1;
 ipg=1;
 
 for j=2:npas,
    y(j)=feval(S,x(:,j))*problem;
    nS=nS+1;

    if y(j) > Ot,
      Ot=y(j);
      ipg=j;
    end
 end

 if ip & n == 2,
   plot([x(1,ipg) x0(1)],[x(2,ipg) x0(2)],'r');
   if ip > 0,
     disp('Pause: hit any key to continue...'); pause;
   end
 end

 xo=p(:,ipg);

 for ig=1:nger,
    rnd1=ones(n,1)*rand(1,npas);
    rnd2=ones(n,1)*rand(1,npas);
    v=w*v+c1*rnd1.*(p-x)+c2*rnd2.*(p(:,ipg)*ones(1,npas)-x);
    x=x+v;
    
    for i=1:n,
       if isb(i),
         j=find(x(i,:) > Ub(i));
         if ~isempty(j),
           x(i,j)=Ub(i);
           v(i,j)=0;
         end
         
         j=find(x(i,:) < Lb(i));
         if ~isempty(j),
           x(i,j)=Lb(i);
           v(i,j)=0;
         end
       end
    end

    for j=1:npas,
       val=feval(S,x(:,j))*problem;
       nS=nS+1;

       if val > y(j),
         y(j)=val;
         p(:,j)=x(:,j);
         
         if val > Ot, 
           Ot=val;
           ipg=j;
         end
       end
    end

    if ip & n == 2,
      plot([x(1,ipg) xo(1)],[x(2,ipg) xo(2)],'r');
      if ip > 0,
        disp('Pause: hit any key to continue...'); pause;
      end
    end

    xo=p(:,ipg);

    w=w-tw;
 end

 if ip & n == 2,
   h2=plot(xo(1),xo(2),'r*');
   legend([h1,h2],'start point','optimum');
 end

 Ot=Ot*problem;
