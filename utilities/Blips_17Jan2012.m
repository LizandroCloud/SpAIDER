function res = Blips(valK,w0,C,L,maxlev)
Pgrad=0;
 for ih=2:(length(w0))-1
    Dgrad(ih-1)= abs(w0(ih))-abs(w0(ih-1));
 end
 il=1;
 normGrad=norm(w0(1:(length(w0))))/(length(w0));
 
 for iw=1:length(Dgrad)
%      Dgsin=Dgrad(iw)*Dgrad(iw-1);
     if abs(Dgrad(iw))>(valK)*normGrad
         Cgrad(il)= iw+1;
         il=il+1;
     end
 end
k=1;
for j=1:1
chg = detcoef(C,L,j);
for i=1:length(chg)
    if(chg(i))~=0
        Pgrad(k)=(i*(2^j)) - 1;
%         Pgrad(k+1)=i*2;
        k=k+1;
    end
end
chg=[]; 
end  


if isempty(Pgrad) 
    
    k=1;
for j=1:2
chg = detcoef(C,L,j);
for i=1:length(chg)
    if(chg(i))~=0
        Pgrad(k)=(i*(2^j)) ;
%         Pgrad(k+1)=i*2;
        k=k+1;
    end
end
      if isempty(Pgrad)
        Pgrad_tr=0;
      else
        Pgrad_tr = unique(Pgrad);  
      end
end
else
Pgrad_tr = unique(Pgrad);
end

if Pgrad==0 
    
    k=1;
for j=1:2
chg = detcoef(C,L,j);
for i=1:length(chg)
    if(chg(i))~=0
        Pgrad(k)=(i*(2^j)) ;
%         Pgrad(k+1)=i*2;
        k=k+1;
    end
end
      if isempty(Pgrad)
        Pgrad_tr=0;
      else
        Pgrad_tr = unique(Pgrad);  
      end
end
else
Pgrad_tr = unique(Pgrad);
end



    res=Pgrad_tr;
 %   res = Cgrad;
end