function [u] = ver_reg_est0(ws,yi,t)

% If the 'est' is 0

for i = 1:length(ws(:,1)) 
       u(i)= ppval(yi,t);     
end

end