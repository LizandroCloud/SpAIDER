
optNLP = optimset( 'Algorithm', 'active-set', 'GradObj', 'off', 'GradConstr', 'off',...
 'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', 10^(-10),...
 'TolFun', 10^(-10), 'TolCon', 10^(-10), 'MaxFunEval', 15000000, 'MaxIter', 5e10);


% options = gaoptimset(@ga) ;


% 'rigrsure' use principle of Stein's Unbiased Risk.
% 'heursure' is an heuristic variant of the first option.
% 'sqtwolog' for universal threshold sqrt(2*log(.)).
% 'minimaxi' for minimax thresholding.
 kh=0; klm=0;
 woptA=[1.51124224601781,0,0,0,0,0,0,0.983623136013518,0,0.883761557525967,0,0.785578905934814,0.732179993265733,0.675441662755597,0.616183501762883,0.482386442306509,0.306097570742583,0.220828520500383,0.207826552552209,0.203979644457102,0.227571093042759,0.251103246197265,0.276487948584357,0,0,0.412676951057678,0.482888472175955,0,0.578319634121656,0.639321912762999,0.702347519033720,0.745418429101277,0.789326840053068,0.834079581125497,0.879422744893521,0.925166780431499,0.971444500295401,1.01827492708295,1.06540568104497,1.11265639458444,1.16017924479922,1.20780213349765,1.26953868355384,1.34586463364483,1.45382063098834,1.58443364180215,1.70192183440081,1.80686624101251,1.90062000682901,1.98427110221374,2.05843044697032,2.12400749090914,2.18230433708268,2.23438250522474,2.28095043755944,2.32283701309356,2.36108371296456,2.39655436674282,2.34172283782656,2.20496143083188,2.07207482108627,1.94293575885235,1.80239736890327,1.65388908164788,1.51833592982289,1.39408188752987,1.28046803130227,1.17658701239669,1.08101829560259,0.992757049058661,0.911432368576592,0.836474067204164,0.766893396825525,0.701980591376041,0.641452854859259,0.584887599206156,0.531517049079258,0.480763274682789,0.432371199954182,0.385973993016210,0.184946878867261];
 woptA = woptA(1:32);
  error0 = randn(1,length(woptA))*0.2;
  ib=1;
%  wopt2 = woptA + error*factor;
%  wopt = wopt2;
 for factor=0.2:0.05:1
     for hgb = 0.2:0.05:0.9
     kh=kh+1; klm=klm+1;
 wname = 'db1';  % tipo de base wavelets.
 tptr = 'cvpshrink'; %cvpshrink
 sorh = 'h';
 scal = 'sln';
 is0 = 2;
 thrs = [];
%  woptA=[1.51124224601781,0,0,0,0,0,0,0.983623136013518,0,0.883761557525967,0,0.785578905934814,0.732179993265733,0.675441662755597,0.616183501762883,0.482386442306509,0.306097570742583,0.220828520500383,0.207826552552209,0.203979644457102,0.227571093042759,0.251103246197265,0.276487948584357,0,0,0.412676951057678,0.482888472175955,0,0.578319634121656,0.639321912762999,0.702347519033720,0.745418429101277,0.789326840053068,0.834079581125497,0.879422744893521,0.925166780431499,0.971444500295401,1.01827492708295,1.06540568104497,1.11265639458444,1.16017924479922,1.20780213349765,1.26953868355384,1.34586463364483,1.45382063098834,1.58443364180215,1.70192183440081,1.80686624101251,1.90062000682901,1.98427110221374,2.05843044697032,2.12400749090914,2.18230433708268,2.23438250522474,2.28095043755944,2.32283701309356,2.36108371296456,2.39655436674282,2.34172283782656,2.20496143083188,2.07207482108627,1.94293575885235,1.80239736890327,1.65388908164788,1.51833592982289,1.39408188752987,1.28046803130227,1.17658701239669,1.08101829560259,0.992757049058661,0.911432368576592,0.836474067204164,0.766893396825525,0.701980591376041,0.641452854859259,0.584887599206156,0.531517049079258,0.480763274682789,0.432371199954182,0.385973993016210,0.184946878867261];
%  woptA = woptA(1:32);
 
ic=1;
 error = error0*factor;
 wopt2 = woptA + error;
 wopt = wopt2;
 maxlev = wmaxlev(length(wopt),wname); 
   j=3;
     [C,L] = wavedec(wopt,maxlev,wname);  % C sao os detalhes...
     [Cw,L] = wavedec(woptA,maxlev,wname); 
     [Cpi,L] = wavedec(error,maxlev,wname);
%    if strcmp(tptr,'cvpshrink')  % bayeshrink
%        [C,L] = wavedec(wopt,maxlev,wname);  % C sao os detalhes...
       i=(maxlev);
        xL = min(abs(C));
%         xL = 0;
        xU = max(abs(C));
%         xU = 2;
        th_objective = @(x)obthr(wopt(ic,:),wname,x);
        m=1;    
        for kh=(xL) :((xU)-(xL) )/99:xU 
%         XX = [xL :(xU-xL )/99:xU];
        [  JJ(m) J1(m) J2(m)] = schaffer(wopt,wname,kh,j,hgb);
%           alpha(m) = Cid(m)^2 / (Cid(m)^2 + (var(error))^2);
           trs(m) = kh;
        m=m+1;
        end
%         plot(trs,JJ); 
%         hold on;
           minJ = min(JJ);
           minJpos = find(JJ==minJ);
           jK = min(minJpos);
           jkmax = max(minJpos);
           minth = trs(jK); 
           maxth = trs(jkmax);
            alfa=ones(1,32)*1e3;
            sum1 = 0; sum2 = 0;
            minJJ(klm) = minJ;
            for  i = L(j+2)+1 : L(j+3)
                Num(j) = Cw(i)*( Cw(i)+Cpi(i)) + sum1;
                sum1 = Num(j);
                Nx(i) = Cw(i)*( Cw(i)+Cpi(i));
                Den(j) = ( ( Cw(i)+Cpi(i))^2 +1 ) + sum2;
                sum2 = Den(j);
                Dx(i) = ( ( Cw(i)+Cpi(i))^2 +1 );
                alfa(i) = Nx(i)/Dx(i);
                thres(i) = ( abs(Cpi(i)));
            end
            
            alpharr(j) = max(thres);
             
            Alpha(klm) =  alpharr(j);
          
 
            TT(klm) = minth;
            TG(klm) = maxth;
            TM(klm) = (maxth + minth) / 2 ;
            
     end    
 
%       for kh=xL :(xU-xL )/99:xU 
% %         XX = [xL :(xU-xL )/99:xU];
%         [  JJ(m) J1(m) J2(m)] = schaffer(wopt,wname,kh,j,hgb);
% %           alpha(m) = Cid(m)^2 / (Cid(m)^2 + (var(error))^2);
%            trs(m) = kh;
%         m=m+1;
%       end
     
     
     
     ktt(ib) = min(abs(TT - ones(1,length(TT))*Alpha(klm)));
     TTB(ib) = Alpha(klm) + ktt(ib);
     ktg(ib) = min(abs(TG - ones(1,length(TG))*Alpha(klm)));
     TGB(ib) = Alpha(klm) - ktg(ib);
     AlphaB(ib) = Alpha(klm);
     
      
     ww = find( (TT)==TTB(ib) );
     if isempty(ww)
        ww = find( (TT)==(Alpha(klm) - ktt(ib) ) );
     end
    
          
    wB(ib) = (min(ww)-1)*0.05 + 0.2;
     
     m=1;
      for kh=abs(xL) :abs(abs(xU)-abs(xL) )/99:xU 
%         XX = [xL :(xU-xL )/99:xU];
        [  JJx(m) J1x(m) J2x(m)] = schaffer(wopt,wname,kh,j,wB(ib));
%           alpha(m) = Cid(m)^2 / (Cid(m)^2 + (var(error))^2);
           trsx(m) = kh;
        m=m+1;
      end
%         figure(2)
%         plot(trsx,JJx); 
%         hold on;
     
     ib = ib+1;
    klm =1;
    e(ib) = factor;
 end     
