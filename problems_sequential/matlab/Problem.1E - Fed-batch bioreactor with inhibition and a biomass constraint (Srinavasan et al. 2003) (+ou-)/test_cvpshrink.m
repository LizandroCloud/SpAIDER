     [C,L] = wavedec(wopt,maxlev,wname);
     h=1; J=0;
     Ci=C;
     
    M=max(  C((L(end-1)+1):end ) ); I=min(  abs(C((L(end-1)+1):end )) );
    for x = I:(abs(M)-abs(I))/10:M
    [C,L] = wavedec(wopt,maxlev,wname);
    tt(h)=x;
    j=1; k=0;

        for i=L(end-1)+1:length(C)
            if abs(C(i))<x;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
                k=j;
                j=j+1;
            end  
        end
    Cth=C;
    w01(1,:) = waverec(Cth,L,wname); % inverse wavelet transformation...
          
       % maximum of detailsss
       j=0; kj=0;
       for i=L(end-1)+1:length(C)
            if abs(C(i))<M;  %fixed threshold criteria
                C(i)=0.0; % forcing these details coefficcients to be zero
                kj=j;
                j=j+1;
            end  
        end
    Cthmax = C;
    w0max(1,:) = waverec(Cthmax,L,wname); % inverse wavelet transformation...
     Jn(h)=norm(wopt - w01)/norm(wopt - w0max);
%     Jn(h)=norm(wopt - w01);
    Jk(h) = (L(end-1)-k)/(L(end-1));
    J(h) = Jn(h) + Jk(h) ;
    h=h+1;
    end
    
    thrsx = ones(1,maxlev)*0.009 ; % threshold value...
    thrs(1,1:length(thrsx),1)=thrsx;
        for i=L(1)+1:length(Ci)
            if abs(Ci(i))<thrsx;  %fixed threshold criteria
                Ci(i)=0.0; % forcing these details coefficcients to be zero
                k=j;
                j=j+1;
            end  
        end
    Cth=Ci;
    w01(1,:) = waverec(Cth,L,wname); % inverse wavelet transformation...
    
    stairs(tt,Jn, 'r'); hold on;   
    stairs(tt,Jk, 'g'); hold on;
    stairs(tt,J, 'b'); hold on;


