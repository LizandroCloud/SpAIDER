function [ u ] = u_teorico(t,umax,sig,ts1,ts2,umin)

eta1 = (1+tanh((t-ts1)/sig))/2;
eta2 = (1+tanh((t-ts2)/sig))/2;
uot = 
u = (1-eta1)*umax + eta1*((1-eta2)*uot+eta2*umin);

end