 
 % Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
 % Algoritmo de Otimização Discretização Total usando Bases Wavelets
 % Adaptativas

 optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 100000);
 close all;
 J = 2; %nível de resolucao...
 x=[];
 d0 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 ];    % Chute inicial 
 dl=[0.0 0.0 0.0 0.0 0.0 0.0 0.0 ];
 du=[1.0 1.0 1.0 1.0 1.0 1.0 1.0 ];
 options = optimset('LargeScale','off');
 tic % contagem do tempo...
%  [d,fval,exitflag,output]  =  fminunc(@myfun,d0,optNLP);
[d, fval, iout,output,lambda,grad,hessian] = fmincon(@(ws)myfun(ws),d0,[],[],[],[],dl,du,[],optNLP);
 toc

 di = d;
 for i=1:length(di) % reorganizando d...
 d(i+1)=di(i);
 end
 d(1)=1/4; % obtido da condição de contorno...
 
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
        
 % construcao da matrix P...
 
 P(1,1,1)=0.5; % definicao...básica
 for n = 0:1:J 
    K = ( 2^(n+1) );
    P(1:K,1:K,K)=[P(1:K/2,1:K/2,K/2)          -1*H(1:K/2,1:K/2,K/2)/(2*K)
                  inv(H(1:K/2,1:K/2,K/2))/(2*K)     zeros(K/2)            ];
 end


% cáculo de x'...
 u=1; %contador para o cálculo...
 for t=0:0.1:2   
 h = haar(t,J);
 x(u) = d*h'; %#ok<AGROW>
 u=u+1;
 end
 
% cáculo de x...
u=1; %contador para o cálculo...

K = ( 2^(J+1) );
coef = d*P(1:K,1:K,K);  % certo...
PP = P(1:K,1:K,K)*P(1:K,1:K,K);  % produto de Hs

%**************************************************
% o=1; 
% for t=0:2/7:2   
%  hg = haar(t,J);
%  dd(o) = t/hg;
%  o=o+1;
%  end

KK = eye(1);
 for n = 0:1:J 
    K = ( 2^(n+1) );
    KK(1:K,1:K)=[KK           zeros(K/2)
                   zeros(K/2)         2/K*eye(K/2)     ];
 end
%**************************************************8

ck=[];
% for i = 1:K
%     ini = 0;
%     for j=1:K
%         ck(i) = P(i,j,K) + ini;
%         ini = ck(i);
%     end
% end
for t=0:.125:1   
h = haar(t,J);
x_opt(u) = d*P(1:K,1:K,K)*h'; %#ok<AGROW>
u=u+1;
tk(u-1)=t;
end


plot(tk,x_opt,'o');


 
 
 
 