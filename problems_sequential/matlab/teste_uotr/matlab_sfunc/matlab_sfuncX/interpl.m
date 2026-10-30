function ureg=interpl(ts,u,ut,h)
idt=max(find(ts>=ut)); idtp=min(find(ts<=ut));
ureg = 0.5*(1+tanh((ts-ut(idt))/h))*u(idt) + 0.5*(1+tanh((ts-ut(idtp))/h))*u(idtp);
%idt=max(find(ts>=ut));
%ureg=u(idt);
end