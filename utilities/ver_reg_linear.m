function [u] = ver_reg_linear(ws,yi,t)

for i = 1:length(ws(:,1)) 

       u(i)= ppval(yi,t);
       
end

end