function [u] = ver_reg(ws,t,te,ks,ns,reg,est)

for i = 1:length(ws(:,1)) 

if est==0
 if reg==1  % regularização ativada
     %  
       u(i) = u_reg1(t,te(i,:),ws(i,:),ks,ns(i));  % regularização ativada

 else
             u(i) = u_reg0(t,te(i,:),ws(i,:),ks,ns(i));  % regularização inativada 
        
           x = te(i,:); 
            y = ws(i,:);
%            yi = interp1(x,[y y(end)],'spline', 'pp');  % interpolacao ativada
%             u(i) = lagrange([y y(end)],x,t);
%             u(i)= ppval(yi,t);
 end

      
else % reg=0
        x = te(i,:); 
        y = ws(i,:);
        yi = interp1(x,[y y(end)],'nearest', 'pp');
        u(i)= ppval(yi,t);
        u(i) = ws(i);
end

end