        function [H M Hi] = H_M(J)
        H=[];
        H(1,1,1)=1;
        for n = J:-1:0 % do maior para o menor...
           K = ( 2^(n+1) );
           for t=1:1:K 
             ti=( t-0.5 )/ K;
             h = haar(ti,n);
            H(1:length(h),t,K)=h;
           end
        end
         Hi = H(1:2^(J+1),1:2^(J+1),2^(J+1));
         % Produto de Hs
         M = H(1:2^(J+1),1:2^(J+1),2^(J+1))*H(1:2^(J+1),1:2^(J+1),2^(J+1));
        end