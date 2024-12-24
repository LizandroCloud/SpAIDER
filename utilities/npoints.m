function [sel_points] = npoints(Pgrad,wopt)
 %*************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrosousa@id.uff.br
 % home page: www.ceopt.com
 % Universidade Federal Fluminense
 % ************************************************************************
 %                          Description
 %
 % Function for evaluating if there exist very closed points.
 %
 % ***********************************************************************
 
 tol = 0.00001;
 k = 1;
 uP = min( length(Pgrad),length(wopt) );
 ind1 = 0; ind2 = 0;
%  Pgrad(1)=0;

 if (length(Pgrad)>2)
 for j=1:uP-1
    difPgrad(j) = abs( wopt(Pgrad(j+1))-wopt(Pgrad(j+1)-1) );
 end
  
 for i=1:length(difPgrad)
     if difPgrad(k) <= tol
         Pgrad(k+1)=0;
         ind1 = 1;
%           if Pgrad(k-1) ~= 0
%             Pgrad(k+1)=k+1;
%             ind1 = 0;
%           else
% %             k = k + 1;  
%           end
     end
%      if k>1
%         if (Pgrad(k-1) ~= 0) && (difPgrad(k)> tol)
%             Pgrad(k+1)= k + 1;
%             ind2 = 1;
%         end
%      end
     
     if k==uP-1
          if ind1 == 1
            Pgrad(k+1)=0;
            break;
          end
     end
     if ind1 ~= 1
         if k<length(difPgrad)
          Pgrad(k+2)= Pgrad(k+1)+1;
          k = k + 2;
          if k>=length(difPgrad)
              break
          end
          
         end
     else
         k = k + 1;
     end
     ind1 = 0; ind2 = 0;
 end
 
 if Pgrad(1) == 1
  Pgrad(1) = 0;
 end
 sel_points = Pgrad(Pgrad~=0);
 else
     sel_points = Pgrad(Pgrad~=0);
 end
end