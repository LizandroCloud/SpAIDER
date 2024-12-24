function [u] = regx(ws,t,te,ks,ns)

for i = 1:length(ws(:,1)) 

       u(i) = u_regS(t,te(i,:),ws(i,:),ks,ns(i));  % regularização ativada

end

end