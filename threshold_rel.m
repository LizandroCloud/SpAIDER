 function [sk,d] = threshold_rel(optSPAIDER,SPAIDER)
 % ************************************************************************
% Dsc. Lizandro de Sousa Santos 
% email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
% home page: www. ????
% Programa de Engenharia Química - COPPE - Federal University of Rio de
% Janeiro -  Brazil
% ************************************************************************
% Description: "threshold_rel" is function to report all calculated thresholds at each iteration.
% inputs: 
% optSPAIDER: inputs of the Spaider algorithm;
% SPAIDER: outputs of Sapaider algorithm;
%
% outputs:
% sk: flag that indicates if everything is ok;
% d: matrix of the reported results.
%
% ************************************************************************
 
ME.message = ''; 
ME2.message = ''; 
d{1,1}=cell2mat({'Results'}); d{1,2}=cell2mat({'of'}); d{1,3}=cell2mat({'Case : '});
d{2,1}=cell2mat({'Control Variables/Threshold Levels'});
d{3,1}=cell2mat({'Iterations'});
        for i=2:length(SPAIDER.nsx(1,1,:))-2;
            d{i+2,1}=i;
        end


jMax=(optSPAIDER.input.nU); % each control variable
w=0;
fi=0;
for j=2:jMax+1 % each control variable
    limit=length(SPAIDER.nsx(j-1,1,:))-3;
    maxTh=length(SPAIDER.thrs(j-1,:,end)) ;
    d{2,w+2}=cell2mat({'Thresholds.u',int2str(j-1)}); % each control variable
        for k=2:maxTh+1
            d{3,w+k}=cell2mat({'level',int2str(k-1)});
        end
    w=j+maxTh-2;    
%         for i=2:length(SPAIDER.thrs(:,j,end));
%             d{i,j}=SPAIDER.nsx(:,j-1,i);
%         end

  
            
   for i=4:limit+3
        intr = fix(((log2(SPAIDER.nsx(j-1,1,i)))));
        ktg = intr;
        for k = 2:ktg + 1
            
            try
                 d{i,k+fi} =  (SPAIDER.thrs(j-1,k-1,i-3) );
            catch ME
            end
            if strcmp(ME.message, 'Index exceeds matrix dimensions.')
                
                
            end
                
         
                 
        end
       
   end
   fi=ktg; 
end

sk=1;
 end
        