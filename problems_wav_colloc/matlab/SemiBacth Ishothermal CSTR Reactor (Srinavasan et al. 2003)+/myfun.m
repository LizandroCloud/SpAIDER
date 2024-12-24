function f = myfun(x) 
% Funcaoobjetivo obtida da colocação wavelets haar...
        f = 1/16 + x(1)^2 + 0.5*x(2)^2 + 0.5*x(3)^2 + (1/4)*x(4)^2 + (1/4)*x(5)^2 + (1/4)*x(6)^2 + ...
           (1/4)*x(7)^2 + (1/8) - x(1)/4 - x(2)/16 - x(3)/16 - x(4)/64 - x(5)/64 - ...
           x(6)/64 - x(7)/64;
       
       J=2;
         % definicacao da matrix haar para todos o níveis...
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
       
       
       
       
         P(1,1,1)=0.5; % definicao...básica
         for n = 0:1:J 
         K = ( 2^(n+1) );
         P(1:K,1:K,K)=[P(1:K/2,1:K/2,K/2)          -1*H(1:K/2,1:K/2,K/2)/(2*K)
                  inv(H(1:K/2,1:K/2,K/2))/(2*K)     zeros(K/2)            ];
         end
         PH(:,:,1)=P(:,:,8)  ;
         K = ( 2^(J+1) );
         PP = P(1:K,1:K,K)*P(1:K,1:K,K);  % produto de Hs
       
         
         KK = eye(1);
            for n = 0:1:J 
                K = ( 2^(n+1) );
                    KK(1:K,1:K)=[KK           zeros(K/2)
                   zeros(K/2)         2/K*eye(K/2)     ];
            end
         
         
         y = [1/4 x];
         f = y*KK*y' + (1/8) - x(1)/4 - x(2)/16 - x(3)/16 - x(4)/64 - x(5)/64 - ...
           x(6)/64 - x(7)/64;
end