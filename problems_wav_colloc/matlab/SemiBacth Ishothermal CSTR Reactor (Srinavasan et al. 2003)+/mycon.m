function [ceq c dceq dc]  = mycon(x,J,d_ini,P,PP,H,Hi,M,K,nU,DWT,iDWT) 

          n=2^(J+1);
%           y(1,1) =99; y(2,1)=99; 
          ko=0;
          for k=1:nU-1
          di = x;
          for i=1+ko:length(x)/nU + ko % reorganizando d...
          y(k,i-ko+1)=x(i);
          end
          ko=i;
          end
          
          for k=3:3
          di = x;
          for i=1+ko:length(x)/nU + ko % reorganizando d...
          y(k,i-ko)=x(i);
          end
          ko=i;
          end
          
          % initial conditions
           y(1,1) = (1 - y(1,2:8)*iDWT(2:8,1))*iDWT(1,1)^-1;
           y(2,1) = (0 - y(2,2:8)*iDWT(2:8,1))*iDWT(1,1)^-1;

            for i=1:3
                xd(i,1:8)=y(i,1:8);
            end    
          
          DX1 = xd(1,:)*Hi;
          DX2 = xd(2,:)*Hi;
          X1 = (xd(1,:))*iDWT + ones(1,8)*0;
          X2 = (xd(2,:))*iDWT ;
          U1 = (xd(3,:))*iDWT;
%           U1 = (x(3,:)*PP)*Hi ;
          g1 = ( ( U1.* (10*X2 - X1 ) ));
          g2 = ((U1.*( (X1 - 10*X2)) - (1-U1).*X2 ));
         
          for i=1:8
          res1(i) = ((DX1(i))/1 - ((g1(i))));
          end
%            res1 = ( res1(1)+4*(res1(2)+res1(4)+res1(6)) + 2*(res1(3)+res1(5)+res1(7)) + res1(8) )*1/3 ;
       
          for i=1:8
          res2(i) = ((DX2(i))/1- ((g2(i))));
          end
           res2 = ( res2(1)+4*(res2(2)+res2(4)+res2(6)) + 2*(res2(3)+res2(5)+ res2(7)) + res2(8) )*1/3 ;
            
          
           c(1:8) = U1(1:8) - 1;
           c(9:16) = -U1(1:8) + 0;
            c(17:24) = X1(1:8) - 1;
            c(25:32) = -X1(1:8) + 0.0000;
           c(33:40) = X2(1:8) - 1;
            c(41:48) = -X2(1:8) + 0.0000;
           c(49:56) = abs(res1) - 0.001;
          c(50) = abs(res2) - 0.001;
          dceq = [];
          dc=[];
          ceq=[];
end