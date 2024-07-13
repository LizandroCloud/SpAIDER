function [sk,d] = spaider_rel(optSPAIDER,SPAIDER)
d{1,1}=cell2mat({'Results'}); d{1,2}=cell2mat({'of'}); d{1,3}=cell2mat({'Case : '});
d{2,1}=cell2mat({'Iterations'});
        for i=2:SPAIDER.iter-1;
            d{i+1,1}=i-1;
        end
jMax=optSPAIDER.input.nU+1;
for j=2:jMax
    d{2,j}=cell2mat({'Stages.u',int2str(j-1)});
        for i=3:SPAIDER.iter;
            d{i,j}=SPAIDER.nsx(j-1,1,i);
        end
end
d{2,jMax+1}=cell2mat({'Acumulated CPU'});
        for i=3:SPAIDER.iter;
            d{i,jMax+1}=SPAIDER.time_ac(i);
        end
        
       
d{2,jMax+2}=cell2mat({'Objective Functional'});
        for i=3:SPAIDER.iter;
            d{i,jMax+2}=abs( SPAIDER.Fobj(i));
        end       
if strcmp(optSPAIDER.wavelets.procedure,'on')        
d{2,jMax+3}=cell2mat({'Convergence Rate'});
        for i=3:SPAIDER.iter-1;
            d{i,jMax+3}=SPAIDER.conv(i);
        end
end
k=1;  

% if strcmp(optSPAIDER.wavelets.procedure,'on')    % MSE
%     for j=jMax+4:jMax-2+jMax+4
%         d{2,j}=cell2mat({'Error.u',int2str(k)});
%             for i=3:SPAIDER.iter-1;  % rows
%                 d{i+1,j}=SPAIDER.mseU(k,i-2) ;
%             end
%             k=k+1;
%     end       
%     else
% end
sk=1;
end
        