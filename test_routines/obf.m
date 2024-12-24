function [ f ] = obf(u_pp1,u_pp0,fi)
    g=0;
    for i=1:length(u_pp1)
        f(i) = (((fi+1)^1)*u_pp0(i) - u_pp1(i))^2 + g;
        g = f(i);
    end
        f = sqrt(g)/length(u_pp1);
       
end