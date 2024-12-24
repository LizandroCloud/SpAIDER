function [Cth, xp, detais1, difd, w01, s, threshold] = wav_thresh_test(C,L,maxlev,wname,detais0,wopt,tptr,sorh,adap,scal,ic,is,is0,normfrac,optSPAIDER)
threshold=ones(1,maxlev);
if adap~=1 % adap =0 (fixed)
    if adap==0
%     Nc = mean(C(L(1)+1:length(C)));
     Nc = optSPAIDER.wavelets.thr
    elseif adap==2
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
else %adp=1     %(adaptive)
     if strcmp(tptr,'gcv')
         [Cth, threshold] = mingcvsoft(C);
         w01(ic,:) = waverec(Cth,L,wname);
         thrsx = ones(1,maxlev)*threshold ; % threshold value...
         thrs(ic,1:length(thrsx),is-is0)=thrsx;
         s=var(C(1:end));
     elseif strcmp(tptr,'bayeshrink')
         
     else
     [w01(ic,:),Cth,L,s,threshold] = wden(wopt(ic,:),tptr,sorh,scal,maxlev,wname); % denoising
     threshold=threshold(end:-1:1);
%      thr = thselect(wopt(ic,:),tptr)*ones(1,maxlev).*s2;  % threshold value...
     thrsx=threshold;  % threshold vector
     end
end

 detais1=Cth;
 difd(ic,1:length(C),is-1)= detais1-detais0;
 
 xp = thrsx;
end