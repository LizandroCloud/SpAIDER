function [u] = ver_reg_est1(ws,yi,t)

    for i = 1:length(ws(:,1)) 
        u(i) = ws(i);
	end

end