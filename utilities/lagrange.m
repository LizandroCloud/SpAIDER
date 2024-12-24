function y = lagrange(wk,xk,x)

    [dum,N] = size(wk);
    y = 0.;
    for k = 1:N
        prod = wk(k);
        for j = 1:N
            if j ~= k
                prod = prod * (x-xk(j)) / (xk(k)-xk(j));
            end
        end
        y = y + prod;
    end
  
end