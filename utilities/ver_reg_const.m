function [u] = ver_reg_const(ws,t,te,ks,ns,reg,est)

for i = 1:length(ws(:,1)) 

       u(i) = u_reg1(t,te(i,:),ws(i,:),ks,ns(i));  % regularização ativada

end

end