function ureg=interpl(ts,u,ut,h)
idt=max(find(ts>=ut)); idtp=min(find(ts<=ut));
if idt>=size(u,1)
    ureg=u(end);
else
    ureg = 0.5*(1+tanh((ts-ut(idt))/h))*u(idt) + 0.5*(1+tanh((ts-ut(idtp))/h))*u(idtp);
end
%idt=max(find(ts>=ut));
%ureg=u(idt);
end