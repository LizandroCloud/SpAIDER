function [Cth, xp, detais1, difd, w01, s, threshold] = wav_thresh(C,L,maxlev,wname,detais0,wopt,tptr,sorh,adap,scal,ic,is,is0,normfrac,optSPAIDER,ns)

optNLP = optimset( 'Algorithm', 'interior-point', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-8),...
 'TolFun', 10^(-8), 'TolCon', 10^(-8), 'MaxFunEval', 15000000, 'MaxIter', 5e10);

threshold=ones(1,maxlev);
if strcmp(adap,'automatic')  % adap =0 (fixed)
     if strcmp(tptr,'gcv') % cross validation
         [Cth, threshold] = mingcvsoft(C);
         w01(ic,:) = waverec(Cth,L,wname);
         thrsx = ones(1,maxlev)*threshold ; % threshold value...
         thrs(ic,1:length(thrsx),is-is0)=thrsx;
         s=var(C(L(2):end));
     elseif strcmp(tptr,'bayes')  % bayeshrink
         sigmahat=median(C)/0.6745;
         threshold=bayes(C(L(2):end),sigmahat);
         Cth = sthresh(C,threshold);
         thrsx = ones(1,maxlev)*threshold ; % threshold value...
         thrs(ic,1:length(thrsx),is-is0)=thrsx;
         w01(ic,:) = waverec(Cth,L,wname);
         s=var(C(L(2):end));
     elseif strcmp(tptr,'cvpshrink')  % bayeshrink
            j=1;
           [C,L] = wavedec(wopt(ic,1:(ns(ic))),maxlev,wname);  % C sao os detalhes...
                   
           for i=(maxlev):-1:min(maxlev,2)
               
                   i=(maxlev);
        xL = min((C));
%         xL = 0;
        xU = max((C));
%         xU = 2;
        th_objective = @(x)obthr(wopt(ic,:),wname,x);
        m=1;
        for kh=xL :(xU-xL )/99:xU 
%         XX = [xL :(xU-xL )/99:xU];
        [  JJ(m) J1(m) J2(m)] = schaffer(wopt(ic,:),wname,kh,j);
%           alpha(m) = Cid(m)^2 / (Cid(m)^2 + (var(error))^2);
           trs(m) = kh;
        m=m+1;
        end
           minJ = min(JJ);
           minJpos = find(JJ==minJ);
           jK = min(minJpos);
           minth = trs(jK);    
               
               
          opt_thres(j) =    minth;
               
               
               
%             [C,L] = wavedec(wopt,maxlev,wname);  % C sao os detalhes...
%             x0 = 0;
%             xL = 0;
%             xU = max(abs(  C((L(end-j)+1):end ) ));
%             th_objective = @(x)obthr(wopt(ic,1:(ns(ic))),wname,x);
% %              [opt_thres, erropt, iout,output,lambda,grad,hessian] = fmincon(  th_objective, x0, [], [], [], [], xL, xU, [] , optNLP);  % opttimal threshold
% %             [opt_thres(j),f(j),exitflag] = ga(@(x)schaffer(wopt(ic,1:(ns(ic))),wname,x,j),1,[],[],[],[],...
% %             xL,xU);
%         
%             [opt_thres(j),f(j),exitflag] = fmincon(@(x)schaffer(wopt(ic,1:(ns(ic))),wname,x,j),x0,[],[],[],[],...
%             xL,xU,[],optNLP);

                         for i=L(end-j)+1:L(end-j+1)
                            if abs(C(i))<opt_thres;  %fixed threshold criteria
                            C(i)=0.0; % forcing these details coefficcients to be zero
                            end  
                         end
                        Cth=C;
                        w01(ic,:) = waverec(Cth,L,wname); % inverse wavelet transformation...
                        s=var(C(1,:));
                  j=j+1;
           end
       thrsx = ones(1,maxlev)*0 ; % threshold value...
%        thrs(ic,1:length(opt_thres),is-is0)=opt_thres;
       thrsx( 1:length(opt_thres) ) = opt_thres;
%        thrsx = thrs;
     else  % sure or visu
     [w01(ic,:),Cth,L,s,threshold] = spaider_wden(wopt(ic,1:(ns(ic))),tptr,sorh,scal,maxlev,wname); % denoising
     threshold=threshold(end:-1:1);
%      thr = thselect(wopt(ic,:),tptr)*ones(1,maxlev).*s2;  % threshold value...
     thrsx=threshold;  % threshold vector
     end
else
     
    if strcmp(adap,'fixed')
%     Nc = mean(C(L(1)+1:length(C)));
     Nc = optSPAIDER.wavelets.fix_threshold
    elseif strcmp(adap,'norm-based')
    Nc = norm(C(L(end-1)+1:end))*normfrac
    end
%     [Cth Nc] = mingcvhard(C(2:end));
    thrsx = ones(1,maxlev)*Nc ; % threshold value...
    thrs(ic,1:length(thrsx),is-is0)=thrsx;
        for i=L(1)+1:length(C)
            if abs(C(i))<Nc;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
            end  
        end
    Cth=C;
    w01(ic,:) = waverec(Cth,L,wname); % inverse wavelet transformation...
    s=var(C(1:end));
 %adp=1     %(adaptive)
    
end

 detais1(ic,1:length(Cth),is-1)=Cth;
 difd(ic,1:length(C),is-1)= detais1(ic,1:length(C),is-1)-detais0(ic,1:length(C),is-1);
 
 xp = thrsx;
end