function f = myfun(x,J,d_ini,P,PP,H,Hi,M,K,nU,DWT,iDWT) 

          

          n=2^(J+1);

            for i=1:nU
                xd(i,1:n)=x(i,1:n);
            end        
          DX1 = xd(1,1:n)*inv(Hi);
          DX2 = (0.5.*(DX1*iDWT).^2)*inv(Hi);
         
          X1 = DX1*iDWT +1;
          X2 = DX2*iDWT ;
        
          f = (X2(end));

         
end