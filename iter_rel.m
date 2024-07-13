function [sk,d] = iter_rel(optSPAIDER,SPAIDER,iter)
% d{1,1}=cell2mat({'Results'}); d{1,2}=cell2mat({'of'}); d{1,3}=cell2mat({'Case : '});
jMax=optSPAIDER.input.nU; k=1;
if strcmp(optSPAIDER.wavelets.procedure,'on')
for j=1:jMax
         d{1,k}=cell2mat({'Details','0'}); d{1,k+1}=cell2mat({'Details','1'});
         for i=2:(SPAIDER.ld0(j,iter))+1;
            d{i,k}=SPAIDER.detais0(j,i-1,iter);
            d{i,k+1}=SPAIDER.detais1(j,i-1,iter);
        end
k=k+2;
end
end
% k=jMax+2;
for j=1:jMax
        d{1,k}=cell2mat({'Time'}); d{1,k+1}=cell2mat({'Control u_',int2str(j)});
        d{1,k+3}=cell2mat({'Time'}); d{1,k+4}=cell2mat({'Thresholded control u_',int2str(j)});
        for i=2:(SPAIDER.nsx(j,:,iter+3))+1;
            d{i,k}=SPAIDER.tx(j,i-1,iter+3);
            d{i,k+1}=SPAIDER.wx(j,i-1,iter+3);
        end
            d{i+1,k}=SPAIDER.tx(j,i,iter+3);
            d{i+1,k+1}=SPAIDER.wx(j,i-1,iter+3);
 
         for i=2:(SPAIDER.nsx(j,:,iter+3))+1;
            d{i,k+3}=SPAIDER.tx(j,i-1,iter+3);
            if strcmp(optSPAIDER.wavelets.procedure,'on') 
            d{i,k+4}=SPAIDER.uth(j,i-1,iter);
            end
         end
        
         d{i+1,k+3}=SPAIDER.tx(j,i,iter+3);
         if strcmp(optSPAIDER.wavelets.procedure,'on')
         d{i+1,k+4}=SPAIDER.uth(j,i-1,iter);
         end
k=k+5; 
end


if strcmp(optSPAIDER.wavelets.procedure,'on')
limit = (SPAIDER.ld0(j,iter))+4;
k=1;
for j=1:jMax
    d{limit,k}=cell2mat({'Mean of detail0: ','u',int2str(j)}); d{limit,k+1}=SPAIDER.md0(j,iter);
    d{limit+1,k}=cell2mat({'Mean of detail1: ','u',int2str(j)}); d{limit+1,k+1}=SPAIDER.md1(j,iter);
    k = k+2;
end   
k=1;
for j=1:jMax
    d{limit+2,k}=cell2mat({'Var of detail0: ','u',int2str(j)}); d{limit+2,k+1}=SPAIDER.md0(j,iter);
    d{limit+3,k}=cell2mat({'Var of detail1: ','u',int2str(j)}); d{limit+3,k+1}=SPAIDER.md1(j,iter);
    k = k+2;
end  

end
sk=1;
end
        