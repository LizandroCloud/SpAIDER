function res = Blips_fix(w0,C,L,maxlev,normfrac)
 
    Normm = norm(C(2:end))*normfrac
    for ih=2:(length(C))
        if abs(C(ih))<Normm
           C(ih) =0.0;
        else
           
        end
    end

Pgrad=0;
%  for ih=2:(length(w0))-1
%     Dgrad(ih-1)= abs(w0(ih))-abs(w0(ih-1));
%  end
%  il=1;
%  normGrad=norm(w0(1:(length(w0))))/(length(w0));
%  
%  for iw=1:length(Dgrad)
% %      Dgsin=Dgrad(iw)*Dgrad(iw-1);
%      if abs(Dgrad(iw))>(valK)*normGrad
%          Cgrad(il)= iw+1;
%          il=il+1;
%      end
%  end
k=0;
for j=1:3  % high frequencies levels
chg = 0;

chg = detcoef(C,L,j);
ko=k+1;

for i=1:length(chg)
%     ko=0;
        gi=1;
    if(chg(i))~=0
        if j==1
        for k=ko:2^j+ko-1
            
            Pgrad(k)=(i*(2^j)) - 2^j + (gi);  % upsampling
%         Pgrad(k+1)=(i*(2^j));
%         k=k+j*2;
            gi=gi+1;
        end
        else
                  for k=ko:(2^j-1):2^j+ko-1
            
            Pgrad(k)=(i*(2^j)) - 2^j + (gi);
%         Pgrad(k+1)=(i*(2^j));
%         k=k+j*2;
            gi=gi+3;
                  end
        end
        ko=k+1;

    end
    
end
chg=[]; 
end  


if isempty(Pgrad) 
    
k=0;
for j=1:3
chg = 0;
chg = detcoef(C,L,j);
ko=k+1;

for i=1:length(chg)
%     ko=0;
gi=1;
    if(chg(i))~=0
        for k=ko:2^j+ko-1
            
            Pgrad(k)=(i*(2^j)) - 2^j + (gi);
%         Pgrad(k+1)=(i*(2^j));
%         k=k+j*2;
            gi=gi+1;

        end
        ko=k+1;
    end
    




end
chg=[]; 
end  % low frequencies levels

else
Pgrad_tr = unique(Pgrad);
end

if Pgrad==0 
    
    k=1;
for j=1:3
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



    res=Pgrad_tr(Pgrad_tr~=0);


end