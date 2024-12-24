function [P PP PH] = P_M(J,H)
    
    

    P(1,1,1)=0.5; % definicao...básica
    for n = 0:1:J 
    K = ( 2^(n+1) );
    P(1:K,1:K,K)=[P(1:K/2,1:K/2,K/2)          -1*H(1:K/2,1:K/2,K/2)/(2*K)
                  inv(H(1:K/2,1:K/2,K/2))/(2*K)     zeros(K/2)            ];
    end
    
     PH(:,:,1)=P(1:2^(J+1),1:2^(J+1),2^(J+1))  ;
     K = ( 2^(J+1) );
     PP = P(1:K,1:K,K)*P(1:K,1:K,K);  % produto de Ps
end