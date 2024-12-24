function f = myfun(x,J,d_ini,P,PP,H,Hi,M,K,nU,DWT,iDWT) 

          
%           y(1,1) = 1000; y(2,1)=8008; y(3,1)=1700;
%           y(1,1) = y(1,2)*iDWT(1,2);
          n=2^(J+1);
%            y(1,1) =99; y(2,1)=99; 
           ko=0;
           for k=1:nU-1
%           di = x;
           for i=1+ko:length(x)/nU + ko % reorganizando d...
          y(k,i-ko+1)=x(i);
          end
%           ko=i;
           end
          
          for k=1:3
          di = x;
          for i=1+ko:length(x)/nU + ko % reorganizando d...
          y(k,i-ko)=x(i);
          end
          ko=i;
          end
          
           y(1,1) = (1 - y(1,2:4)*iDWT(2:4,1))*iDWT(1,1)^-1;
            y(2,1) = (0 - y(2,2:4)*iDWT(2:4,1))*iDWT(1,1)^-1;
%             y(3,1) = (1 - y(3,2:4)*iDWT(2:4,1))*iDWT(1,1)^-1;
%            x=y;
          %
% n=2^(J+1);
%           ko=0;
%           for k=1:nU
%           di = x;
%           for i=1+ko:length(x)/nU + ko % reorganizando d...
%           y(k,i-ko)=x(i);
%           end
% %          d_ini(k) = -y(k,2)*Hi(2,1)*PP(2,1) - y(k,3)*Hi(3,1)*PP(3,1) - y(k,5)*Hi(5,1)*PP(5,1);
% %          d_ini(k,2) = -y(k,2)*Hi(2,1)*PP(2,1) - y(k,3)*Hi(3,1)*PP(3,1) - y(k,5)*Hi(5,1)*PP(5,1);
% %          y(k,1)=d_ini(k); % obtido da condição de contorno...
%           ko=i;
%           end
            for i=1:3
                xd(i,1:4)=y(i,1:4);
            end        
          DX1 = xd(1,1:4)*Hi;
          DX2 = xd(2,1:4)*Hi;
           X1 = (xd(1,:))*iDWT + ones(1,8)*0;
          X2 = (xd(2,1:4))*iDWT ;
          U1 = (xd(3,1:4))*iDWT;
%           U1 = (x(3,:)*PP)*Hi ;
          g1 = ( ( U1.* (10.*X2 - X1 ) ));
          g2 = ((U1.*( (X1 - 10.*X2)) - (1-U1).*X2 ));
         res1 = (DX1(end)+DX1(1))/4 - ((g1(end))-1)/1;
         res2 = (DX2(end)+DX2(1))/8 - ((g2(end))-0)/1;
          sum(DX1)*0.1;
          ix1 = (g1(end)-1)*0.5;
          ix2 = (g2(end)-0)*0.5;
          f = -1+ ( X1(end))+ (X2(end));
%           f = -1+(y(1,:))*iDWT + (y(2,:))*iDWT ;
          f = f(end);
end