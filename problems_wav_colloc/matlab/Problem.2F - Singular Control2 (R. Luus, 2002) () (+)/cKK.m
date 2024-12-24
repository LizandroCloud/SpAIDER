function [KK] = cKK(J)

% Integral do produto de Hs
         KK = eye(1);  % integral de M
            for n = 0:1:J 
                K = ( 2^(n+1) );
                    KK(1:K,1:K)=[KK           zeros(K/2)
                   zeros(K/2)         2/K*eye(K/2)     ];
            end
            
end