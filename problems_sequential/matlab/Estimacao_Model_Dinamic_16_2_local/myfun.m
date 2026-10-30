function f = myfun(v,modelo)

H=diag(2*ones(length(modelo.c),1));

f=1/2*(v'*H*v);