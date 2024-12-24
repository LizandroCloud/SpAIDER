% function J = schaffer(wopt,wname,x,jo)
 
%  if nargout == 1   % only objective functional...
    wopt=[ 1 1 1 1 3 3 3 3 3 3 3 3 -1 -1 -1 -1 -1 -1 6 6 6 6 6 6 1 1 1 1 1 7 7 7];
    x= randn(1,32);
    woptx = wopt + x*0.1;
    wname='db1';
    jo=1;
    j=1; k=0;
    [cfd,detais0,maxlev,C,L] = wav_mesh(wopt,3,3,wname);
     M=max(  C((L(end-jo)+1):end ) ); I=min(  C((L(end-jo)+1):end ) );
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
    J(1)=norm(wopt - w01)/norm(wopt - w0max);
    J(2) = ((L(end)-jo)-k)/(L(end)-jo);
    J = J(1) + J(2) ;
    
%      else % objective function and constraints...
%  end

% end