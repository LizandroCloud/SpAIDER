function res = Blips(valK,w0,C,L)

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
 
 % pontos prospectivos...
 %Pgrad(1)= Cgrad(1)+1;
%   Pgrad(2)= Cgrad(1)+2;
 %Pgrad(2)= Cgrad(length(Cgrad))-2;
%   Pgrad(4)= Cgrad(2);
 Pgrad=Cgrad;
 res=Pgrad;
end