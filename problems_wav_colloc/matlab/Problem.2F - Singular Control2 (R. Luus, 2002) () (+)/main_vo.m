 
 % Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
 % Algoritmo de Otimização Discretização Total usando Bases Wavelets
 % Adaptativas
%  resetall;
 optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-6), 'TolCon', 10^(-10), 'MaxFunEval', 100000);
 close all;
 
%************************************************** 
wname = 'db1';  % wavelets type
tptr = 'sqtwolog'; % thresholding procedure
sorh = 'h';  % hard or soft threshold option
scal = 'mln'; % scale procedure
Wav = 1;  % 1 - Wav; 2 -> ADE...
adap=1; % adaption on(thresholding wavelets) or off (fixed)
options = optimset('LargeScale','off'); 

%************************************************** 
 J=3; nUt=1;  nU=nUt; ns = (2^(J+1));
 u_ub=1;
 u_lb=0;
 d0=[]; dl=[]; du=[];
%  problem_Bs   % getting states...
%  for i=1:SPAIDER.n_s
%     spline([1:length(SPAIDER.fst(:,i))],d_ini(:,i),1:8)
%  end
%  
%   for i=1:nu
%     spline([1:length(SPAIDER.wopt(:,i))],d_ini(:,i+SPAIDER.n_s),1:8)
%  end
 d_ini(1:16,1)=[-1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 1];
 d_ini(1:16,2)=[0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16 0.16];
 d_ini(1:16,3)=[1 1 1 1 0.5 0.5 0.5 1 1 1 1 0.0 0.0 0.0 0 0];
      for i=1:length (    d_ini(1,:) )
%      d_inis(1:16,i)= spline(1:length(    d_ini(:,i) ),d_ini(:,i),1:16)
%       Y(:,i) = interp(d_ini(:,i),0)
     X(:,i) = upsample( d_ini(:,i),1 );
      end
%  
 d_ini=X;
 u_lb = d_ini - 0.1*d_ini ;
 u_ub = d_ini + 0.1*d_ini ;
 u0 = d_ini;
%  u_lb(:,3) = ones(ns,1)*-10;
%  u_ub(:,3) = ones(ns,1)*10;
  u_lb(:,1) = ones(ns,1)*-1;
  u_ub(:,1) = ones(ns,1)*1;
%  u_lb(:,2) = ones(ns,1)*-10;
%  u_ub(:,2) = ones(ns,1)*10;
%  J = 3;



 ns(1:3) = (2^(J+1)); % initial number of stages
 x=[];
%  u_lb=[0.4 0.4 0]; u_ub=[0.6 0.6 1];  u0= [0.5 0.5 0.5];
[H M Hi] = H_M(J);
[P PP PH] = P_M(J,H);
[KK] = cKK(J);
[DWT iDWT] = wav_mat(PH, Hi ,J);
 for ic = 1:nUt 
     d_ini(:,ic)=u0(:,ic)/(PP(:,1)'*Hi(:,1));
 end
 for ji = 1:nUt; % for each control, variable
        ns(ji) = (2^(J+1)); % initial number of stages
        dL(ji,1:ns(ji)) = u_lb(:,ji);  % lower boundary of u
%         dL(ji,1:8) = u_lb(:,ji)'*inv(Hi)*inv(eye(8,8))*inv(PP)*inv(eye(8,8));  
        dU(ji,1:ns(ji)) = u_ub(:,ji);  % upper boundary of u
        d0(ji,1:ns(ji)) = u0(:,ji);  % initial estimative da u...
%         d0(ji,1:ns(ji)) = ( (u0(ji)*ones(ns(ji),1))'*inv(Hi))*inv(PP);
%         d0(ji,1:ns(ji)) = ((PP\( u0(ji)*ones(ns(ji),1) ))')/Hi;
        axc = dL(ji,:)*DWT;
%         ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
%         w_frozen = w0;
%         check = [2 3 4 5 6 7 8 9];
%         ts_frozen = ts;
%         ts_a = ts;
 end
 
  for ji = 1:nUt; % for each control, variable
        ns(ji) = (2^(J+1)); % initial number of stages
        dL(ji,1:ns(ji)) = u_lb(:,ji);  % lower boundary of u
%         dL(ji,1:8) = u_lb(:,ji)'*inv(Hi)*inv(eye(8,8))*inv(PP)*inv(eye(8,8));  
        dU(ji,1:ns(ji)) = u_ub(:,ji);  % upper boundary of u
        d0(ji,1:ns(ji)) = u0(:,ji);  % initial estimative da u...
%         d0(ji,1:ns(ji)) = ( (u0(ji)*ones(ns(ji),1))'*inv(Hi))*inv(PP);
%         d0(ji,1:ns(ji)) = ((PP\( u0(ji)*ones(ns(ji),1) ))')/Hi;
        axc = dL(ji,:);
%         ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
%         w_frozen = w0;
%         check = [2 3 4 5 6 7 8 9];
%         ts_frozen = ts;
%         ts_a = ts;
 end
 
 
 
%   usup = [ dU(3,:)]/(PP)'*Hi;
%   uinf = [ dL(3,:)]/(PP)'*Hi;
%   uini = [ d0(3,:)]/(PP)'*Hi;
 for ji = 1:nUt; % for each control, variable
        for k=1:ns(ji)
        usup1(k) = max(dU(ji,k),dL(ji,k)) ;
        usup2(k) = min(dU(ji,k),dL(ji,k)) ;
        mA= usup1(k)+0.8*abs(usup1(k));
        mB= usup2(k)-0.8*abs(usup2(k));
        
        if (mA-mB)<0.001
            mA = mA+10.1;
            mB = mB-10.1;
        end
        
        dL(ji,k) = mB;  % lower boundary of u
        dU(ji,k) = mA;  % upper boundary of u
%         d0(ji,1:ns(ji)-1) = u0(ji)*ones(ns(ji)-1,1);  % initial estimative da u...
%         ts(j,1:ns(j)+1) = (t0:(tf-t0)/ns(j):tf); % initial time discretization...
%         w_frozen = w0;
%         check = [2 3 4 5 6 7 8 9];
%         ts_frozen = ts;
%         ts_a = ts;
        end
 end
% % for i=1:3
%  dU(i,:) = usup1(2:end);
%  dL(i,:) = usup2(2:end);
%  end
%  
 
 

 
 k=0; ki=1;
    for i=1:nUt
     for j=k+1:ns(i)+k
         d00n(j)=(d0(i,ki)-dL(i,ki))/(dU(i,ki)-dL(i,ki));
        d00(j)=u0(ki,i);
        d00a(j) = d0(i,ki);
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
   k=0; ki=1;
   for i=1:nUt
     for j=k+1:ns(i)+k
         dLL(j)=u_lb(ki,i);
         dUU(j)=u_ub(ki,i)';
          dLLn(j)=0;
          dUUn(j)=1;
         ki=ki+1;
     end
     k = j;
     ki=1;
   end
%    koo=k; kii=ki;
% %    k=0; ki=1;
%     for i=3:3
%      for j=k+1:ns(i)+k
%          d00n(j)=(d0(i,ki)-dL(i,ki))/(dU(i,ki)-dL(i,ki));
%         d00(j)=u0(ki,i);
%         d00a(j) = d0(i,ki);
%          ki=ki+1;
%      end
%      k = j;
%      ki=1;
%    end
% %    k=0; ki=1;
%    for i=3:3
%      for j=koo+1:ns(i)+koo
%          dLL(j)=u_lb(kii,i);
%          dUU(j)=u_ub(kii,i)';
%           dLLn(j)=0;
%           dUUn(j)=1;
%          kii=kii+1;
%      end
%      koo = j;
%      kii=1;
%    end
 
%  for ic=1:3
%   d0(ic,1:2^(J+1)-1) = 0*ones(2^(J+1)-1,1)';    % Chute inicial 
%   dl(ic,1:2^(J+1)-1)= u_lb(ic)*ones(2^(J+1)-1,1)'; 
%   du(ic,1:2^(J+1)-1)= u_ub(ic)*ones(2^(J+1)-1,1)'; 
%  end
%  options = optimset('LargeScale','off');
 %**************************************************
%  d_ini=[0  0 1/4];
 
 tic
 for j = J:J; %nível de resolucao...
 is=J;
  % contagem do tempo...
%  [d,fval,exitflag,output]  =  fminunc(@myfun,d0,optNLP);
% [H M Hi] = H_M(J);
% [P PP] = P_M(J,H);
% [KK] = KK(J);
% usup = [1 dU(3,:)]*inv(Hi)*inv(PP);

% d_ini(1) = -d0(1,2)*Hi(2,1)*PP(2,1) - d0(1,3)*Hi(3,1)*PP(3,1) - d0(1,5)*Hi(5,1)*PP(5,1);
% d_ini(2) = -d0(2,2)*Hi(2,1)*PP(2,1) - d0(2,3)*Hi(3,1)*PP(3,1) - d0(2,5)*Hi(5,1)*PP(5,1);
[d, fval, iout,output,lambda,grad,hessian] = fmincon(@(ws)myfun_vo(ws,J,d_ini,P,PP,H,Hi,M,KK,nU,DWT,iDWT),d00,[],[],[],[],dLL,dUU,@(ws)mycon(ws,J,d_ini,P,PP,H,Hi,M,KK,nU,DWT,iDWT),optNLP);
 
%**************************************************

  
             for i=1:nU
                xd(i,1:ns)=d(i,1:ns);
            end 
          DX1 = xd(1,1:ns)*inv(Hi);
          DX2 = (0.5*DX1*iDWT.*DX1*iDWT)*inv(Hi);
         
          X1 = DX1*iDWT ;
          X2 = DX2*iDWT ;
%           for k=3:3
%           di = x;
%           for i=1+ko:length(d)/nU + ko % reorganizando d...
%           y(k,i-ko)=d(i);
%           end
%           ko=i;
%           end

%           y(1,1) = (1 - y(1,2:8)*iDWT(2:8,1))*iDWT(1,1)^-1;
%            y(2,1) = (0 - y(2,2:8)*iDWT(2:8,1))*iDWT(1,1)^-1;
%            y(3,1) = (1 - y(3,2:8)*iDWT(2:8,1))*iDWT(1,1)^-1;
          
%              for i=1:3
%                 xd(i,1:4)=y(i,1:4)*iDWT;
%              end

          
 end
 

 
 
toc
% stairs(tk,x_opt,'-');


 
 
 
 