xref = SPAIDER.state;
tref = SPAIDER.time_state;
tplot2 = tplot + 1e-8*rand(length(tplot),1);
tplot3 = sort(tplot2);
tplot3(1) = 0; tplot3(end)=6;
xx = interp1(tplot3,xplot,[0:0.1:6]);
xxp = interp1(tplot3,xplot,'cubic','pp');
xxref = interp1(tref,xref,'cubic','pp');
SS0=0;

% Point Point SSE
for j=1:3
    for i=1:length(xx(:,1))
        SSE = ( xx(i,j) - xref(i,j) )^2;
        SS0 = SS0 + SSE;
    end
end

SSE = SS0;

% Collocation Points SSE

SS0=0;
times = roundn (optimout.t,-3);
times = unique(times);
% for j=1:3
    for i=1:length(times)
        SSE = ( ppval(xxp,times(i)) - ppval(xxref,times(i))  ).^2;
        SS0 = SS0 + SSE;
    end
% end
    
SSE = sum(SS0);


