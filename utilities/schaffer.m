function [J J1 J2]= schaffer(wopt,wname,x,jo)

%  if nargout == 3   % only objective functional...
%      jo=5;

    j=1; k=0;
    [cfd,detais0,maxlev,C,L] = wav_mesh(wopt,3,3,wname);
%     jo = maxlev;
    C_or = C;
     M=max(  C((L(end-jo)+1):end ) ); I=min(  C((L(end-jo)+1):end ) );
    thrsx = ones(1,maxlev)*x ; % threshold value...
    thrs(1,1:length(thrsx),1)=thrsx;
%         for i=L(jo+1)+1:L(jo+2)+L(jo+1)-2
            for i=L(end-jo)+1:L(end-jo+1)
            if abs(abs(C(i)))<x;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
                k=j;
                j=j+1;
            end  
        end
    Cth=C;
    w01(1,:) = waverec(Cth,L,wname); % inverse wavelet transformation...
    s=var(C(1:end));
       
       % maximum of detailsss
       j=0; kj=0;
       for i=L(end-jo)+1:L(end-jo+1)
            if abs(C_or(i))<M;  %fixed threshold criteria
                C_or(i)=0.0; % forcing these details coefficcients to be zero
                kj=j;
                j=j+1;
            end  
        end
    Cthmax = C_or; % unthouched
    w0max(1,:) = waverec(Cthmax,L,wname); % inverse wavelet transformation...
    J(1)=norm(wopt - w01)/norm(wopt - w0max);
    J(2) = ((L(end)-jo)-k)/(L(end-jo));
    J1=J(1); J2=J(2);
%     J = hgb*J1 +  (1-hgb)*J2 ;
    J = J1 +  J2 ;
%      else % objective function and constraints...
%  end

end