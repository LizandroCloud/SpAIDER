function f = myfun(x,J,d_ini,DWT,iDWT,PH) 
% Funcaoobjetivo obtida da colocação wavelets haar...
%          f = 1/16 + x(1)^2 + 0.5*x(2)^2 + 0.5*x(3)^2 + (1/4)*x(4)^2 + (1/4)*x(5)^2 + (1/4)*x(6)^2 + ...
%            (1/4)*x(7)^2 + (1/8) - x(1)/4 - x(2)/16 - x(3)/16 - x(4)/64 - x(5)/64 - ...
%            x(6)/64 - x(7)/64;
       
          di = x;
          for i=1:length(x) % reorganizando d...
          y(i+1)=x(i);
          end
          y(1)=d_ini; % obtido da condição de contorno...
          
       % definicacao da matrix haar para todos o níveis...
        k=1;
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
        K = ( 2^(J+1) );
       for t=1:1:K 
            ti=( t-0.5 )/ K;
            tid(k)=ti;
            k=k+1;
       end
      
%          [cfd,detais0,maxlev,C,L] = wav_mesh(tid,3,1,'db1');
     
         % Matriz operacional P (Haar)
         P(1,1,1)=0.5; % definicao...básica
         for n = 0:1:J 
         K = ( 2^(n+1) );
         P(1:K,1:K,K)=[P(1:K/2,1:K/2,K/2)          -1*H(1:K/2,1:K/2,K/2)/(2*K)
                  inv(H(1:K/2,1:K/2,K/2))/(2*K)     zeros(K/2)            ];
         end
         
         
         % Produto de Ps
         PH(:,:,1)=P(1:2^(J+1),1:2^(J+1),2^(J+1))  ;
         K = ( 2^(J+1) );
         PP = P(1:K,1:K,K)*P(1:K,1:K,K);  % produto de Ps
       
         % Produto de Hs
         M = H(1:2^(J+1),1:2^(J+1),2^(J+1))*H(1:2^(J+1),1:2^(J+1),2^(J+1));
         
         % Integral do produto de Hs
         KK = eye(1);  % integral de M
            for n = 0:1:J 
                K = ( 2^(n+1) );
                    KK(1:K,1:K)=[KK           zeros(K/2)
                   zeros(K/2)         2/K*eye(K/2)     ];
            end
         [H M Hi] = H_M(J);
         [P PP PH] = P_M(J,H);
         [KK] = cKK(J);
         [DWT iDWT] = wav_mat(PH, Hi ,J);
%          A=C*KK;
     
                a1=y*inv(Hi);
                a2=y*iDWT;
         
          f = y*KK*y' + (1/8) - y(2)/4 - y(3)/16 - y(4)/16 - y(5)/64 - y(6)/64 - ...
             y(7)/64 - y(8)/64;
%           f = y*KK*y' ;
end