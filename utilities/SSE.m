xref = SPAIDER.state;
tplot2 = tplot + 1e-15*rand(204,1);
tplot3 = sort(tplot2);
tplot3(1) = 0;
xx=interp1(tplot3,xplot,[0:0.1:6]);
SS0=0;
for j=1:3
    for i=1:length(xx(:,1))
        SSE = ( xx(i,j) - xref(i,j) )^2;
        SS0 = SSE0 + SSE;
    end
end