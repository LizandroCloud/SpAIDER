function [J,dJ] = obthr(wopt,wname,x)

 if nargout == 1   % only objective functional...

    j=1; k=0;
    [cfd,detais0,maxlev,C,L] = wav_mesh(wopt,3,3,wname);
     M=max(  C((L(end-1)+1):end ) ); I=min(  C((L(end-1)+1):end ) );
    thrsx = ones(1,maxlev)*x ; % threshold value...
    thrs(1,1:length(thrsx),1)=thrsx;
        for i=L(1)+1:length(C)
            if abs(C(i))<x;  %fixed threshold criteria
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
       for i=L(1)+1:length(C)
            if abs(C(i))<M;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
                kj=j;
                j=j+1;
            end  
        end
    Cthmax = C;
    w0max(1,:) = waverec(Cthmax,L,wname); % inverse wavelet transformation...
    Jn=norm(wopt - w01)/norm(wopt - w0max);
    Jk = ((L(end)-2)-k)/(L(end)-2);
    J = Jn + Jk ;
    
     else % objective function and constraints...
 end

end