function [ts_ref w_ref ns_ref point] = meshref(Pgrad,ns,ic,ts_a,w02,t0,tf)

% ************************************************************************
 % Dsc. Lizandro de Sousa Santos 
 % email: lizandrossantos@gmail.com or lizandro@peq.coppe.ufrj.br
 % home page: www. ????
 % Programa de Engenharia Qu�mica - COPPE - Federal University of Rio de
 % Janeiro -  Brazil
 % ************************************************************************
 % Function: meshref.m                           
 %                          Description
 % This routine is used to refine the control profile grid. The vector
 % Pgrad indicates the prospective points that will be refined. The vector
 % w02 is the optimal control profile that will be refined. The refinement procedure
 % keeps the the selected point and adds another one in its neighbhood, more
 % specifically in the stage before, in a distance from the select point
 % that will dependend on the weight parameter dref [0 < dref < 1].
 %
 % for more information:
 % http://www.sciencedirect.com/science/article/pii/B9780444595201501032
 % ************************************************************************
 %
 % [ts_ref w_ref ns_ref point] = meshref(Pgrad,ns,ic,ts_a,w02,t0,tf)
 %
 % input variables:
 % Pgrad -> vector of indexes of prospective points: 
 %          i.e. u = [ 0.11 0.2 0.3 1]  
 %          if 0.2 and 0.3 are selected, then Pgrad = [2 3]
 % ns -> number of stages
 % ic -> control variable ic (ic<=1<=nU)
 % ts_a -> vector of time stages
 % w02 -> control profile (grid)
 % t0-> initial time
 % tf -> final time
 %
 % output variables:
 %
 % ts_ref -> vector of refined time stages
 % w_ref -> vector of refined control profile
 % ns_ref -> number of stages of refined mesh
 % point -> position of Pgrad in the new grid
 %
 % ************************************************************************
 % last changes:
 % 10-30-2012: creation...
 

jj=1;
for j=1:1:ns(ic)  % It will be improved in future...
    ind = 0;
        for hs=1:length(Pgrad)  % veryfiing Pgrad...
                delta(hs) = j - Pgrad(hs);
                if delta(hs)==0 
                    ind=1;
                end
                
                if j==1
                 delta2 = j+1 - Pgrad(hs);
                 if delta2 ==0 
                     ind=1;
                 end
                end 
                 if j==ns(ic)
                 delta2 = j+1 - Pgrad(hs);
                 if delta2 ==0 
                     ind=1;
                 end
                 end
                 
        end
                if (ind==1)
                     point(jj)=1;
                     point(jj+1)=1;
                     
                     if j==1
                     ts_ref(ic,jj+1)=ts_a(ic,j+1)-(ts_a(ic,j+1)- 0)*0.4; 
                     ts_ref(ic,jj+2)=ts_a(ic,j+1);
                     w_ref(ic,jj)=w02(ic,j);  w_ref(ic,jj+1)=w02(ic,j);      %mesh refinement
                     elseif j==ns(ic)
                     ts_ref(ic,jj+1)=ts_a(ic,j+1)-(ts_a(ic,j+1)-ts_a(ic,j))*0.4; 
                     ts_ref(ic,jj+2)=ts_a(ic,ns(ic)+1);
                     w_ref(ic,jj)=w02(ic,j);  w_ref(ic,jj+1)=w02(ic,j);    %mesh refinement
                     else
                     ts_ref(ic,jj+1)=ts_a(ic,j+1)-(ts_a(ic,j+1)-ts_a(ic,j))*0.4; 
                     ts_ref(ic,jj+2)=ts_a(ic,j+1);
                     w_ref(ic,jj)=w02(ic,j);  w_ref(ic,jj+1)=w02(ic,j);     %mesh refinement
                     end
                     jj=jj+2;
                else
                     w_ref(ic,jj)=w02(ic,j);      %mesh refinement
                     ts_ref(ic,jj+1)=ts_a(ic,j+1);
                     point(jj)=0;
                     point(jj+1)=0;
                     jj=jj+1;
                end   
end   
tam = 0;
 for i=1:length(ts_ref(ic,:))  % constructing the new grid (multivariable control)
     diff = tf-ts_ref(ic,i);
     if (diff) <= 0.001
         tam = i;
         ns_ref(ic) = tam-1;
         break;
     end 
 ns(ic) = 2;
 
 end
end