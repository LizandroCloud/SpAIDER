 
 % Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
 % Algoritmo de Otimização Discretização Total usando Bases Wavelets
 % Adaptativas

 optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 100000);
 close all;
 
%************************************************** 
wname = 'db1';  % wavelets type
tptr = 'sqtwolog'; % thresholding procedure
sorh = 'h';  % hard or soft threshold option
scal = 'mln'; % scale procedure
Wav = 1;  % 1 - Wav; 2 -> ADE...
adap='fixed'; % adaption on(thresholding wavelets) or off (fixed)
 
%************************************************** 
 nUt=1
 u_ub=1;
 u_lb=0;
 d0=[]; dl=[]; du=[];
 for ic = 1:nUt 
 end   
 J = 2;
 
 [H M Hi] = H_M(J);
[P PH] = P_M(J,H);
[KK] = cKK(J);
[DWT iDWT] = wav_mat(PH, Hi ,J);
%  for ic = 1:nUt 
%      d_ini(ic)=u0(ic)/(PH(:,1)'*Hi(:,1));
%  end
 
 
 x=[];
 d_ini=1/4;
 d0(ic,1:2^(J+1)-1) = 1*ones(2^(J+1)-1,1)';    % Chute inicial 
 dl(ic,1:2^(J+1)-1)= u_lb*ones(2^(J+1)-1,1)'; 
 du(ic,1:2^(J+1)-1)= u_ub*ones(2^(J+1)-1,1)'; 
 options = optimset('LargeScale','off');
 tic
 for J = 2:2; %nível de resolucao...
 is=J;
  % contagem do tempo...
%  [d,fval,exitflag,output]  =  fminunc(@myfun,d0,optNLP);
[d, fval, iout,output,lambda,grad,hessian] = fmincon(@(ws)myfun_vo(ws,J,d_ini,DWT,iDWT,PH),d0,[],[],[],[],dl,du,[],optNLP);
 

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


KK = eye(1);
 for n = 0:1:J 
    K = ( 2^(n+1) );
    KK(1:K,1:K)=[KK           zeros(K/2)
                   zeros(K/2)         2/K*eye(K/2)     ];
 end
%**************************************************


 maxlev=J+1; is=4; is0=3;
 [cfd,detais0,maxlev,C,L] = wav_mesh(d,is,1,'db1');
 [Cth, thrs1, detais1, difd, w01, s, threshold] = wav_thresh(d,L,maxlev,wname,detais0,d,tptr,sorh,adap,scal,ic,is,2,0.5); % Wavelets thresholding
            thrs(ic,1:length(thrs1),is-is0)=thrs1; 
%             [Pgrad] = adapt(w01,Cth,L,maxlev,ic,tunn,ts,adap,optSPAIDER.wavelets.normfrac) % prospective points...
       
%**************************************************

Pgrad=2:length(d(: ));
ns=2^(J+1);
ic=1;
ts_a=[0:1/(ns):1]; ts=ts_a;
w02=d;
wfact=0.5;
t0=0;
tf=1;
is=3;
 [Pgradx] = adapt(w01,Cth,L,maxlev,ic,32,ts,adap,0.5) % prospective points...

 [ts_refu w_refu ns_refu point] = meshref(Pgradx,ns,ic,ts_a,w02,t0,tf,wfact);
 [ts_ref w_ref ns_ref point] = meshref(Pgrad,ns,ic,ts_a,w02,t0,tf,wfact);
 
 ts(ic,1:ns_ref+1)=ts_ref;
 w0ns(ic,1:ns_ref)=w_ref;
 ns(ic)=ns_ref;  

%**************************************************
du(ic,1:ns(ic)-1) = u_ub(ic)*ones(ns(ic)-1,1);  % lower boundary of u
dl(ic,1:ns(ic)-1) = u_lb(ic)*ones(ns(ic)-1,1);  % upper boundary of u
d0(ic,1:ns(ic)-1)  =w_ref(ic,2:( ns (ic) ) );

 end

%**************************************************

for t=0:.125/J:1   
h = haar(t,J);
x_opt(u) = d*P(1:K,1:K,K)*h'; %#ok<AGROW>
u=u+1;
tk(u-1)=t;
end

toc
stairs(tk,x_opt,'-');


 
 
 
 