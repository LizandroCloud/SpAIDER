function  [PgradN dtA] = Pgrad_ref(Pgrad,spp,ts,ic)
Pgrad=Pgrad(Pgrad~=0);
if length(Pgrad)==1
    i=1;
    dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % control cosntraint
    if ( dtA(i)) >ts(ic,end)/spp
    PgradN(i) = Pgrad(i);
    else  % evaluating neighbors points
%         Pgrad(i) = Pgrad(i)+1;
%         dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % neq dtA
%             if ( dtA(i)) > ts(ic,end)/spp
%                 PgradN(i) = Pgrad(i);
%                     else  % evaluating neighbors points
%                          Pgrad(i) = Pgrad(i)-2;
%                                 dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % neq dtA
%                                     if ( dtA(i)) > ts(ic,end)/spp
%                                         PgradN(i) = Pgrad(i);
%                                     else
                                        PgradN(i) = 0;
    end
else
 biaso=0; bias=0;
 for j=1:length(Pgrad)   
    if Pgrad(j)>length(ts(ic,:))
      bias=1 + biaso;
    end
      biaso=bias;
 end
    bias
 for i=1:(length(Pgrad)-bias)-1
% for i=1:length(ts(ic,:))-1
    dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % control cosntraint
    if ( dtA(i)) >ts(ic,end)/spp
    PgradN(i) = Pgrad(i);
        if i==length(Pgrad)-1
             PgradN(i+1) = Pgrad(i+1);
        end
    else  % evaluating neighbors points
%         Pgrad(i) = Pgrad(i)+1;
%         dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % neq dtA
%             if ( dtA(i)) > ts(ic,end)/spp
%                 PgradN(i) = Pgrad(i);
%                     else  % evaluating neighbors points
%                          Pgrad(i) = Pgrad(i)-2;
%                                 dtA(i) = ts(ic,Pgrad(i)+1) - ts(ic,( Pgrad(i)) );  % neq dtA
%                                     if ( dtA(i)) > ts(ic,end)/spp
%                                         PgradN(i) = Pgrad(i);
%                                     else
                                        PgradN(i) = 0;
                                    end
            
end   
 PgradN = unique(PgradN);  
end