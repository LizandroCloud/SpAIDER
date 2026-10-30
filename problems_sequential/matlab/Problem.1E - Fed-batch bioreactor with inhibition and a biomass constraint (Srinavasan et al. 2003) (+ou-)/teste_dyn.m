
u=1;
% x = [1 2 2 3];
[ dx par ] = dyn_model2;
parameters=char(par);
eq = char(dx);  % entraria com parameters e eq e faria o eval dentro...

% syms x(1) x(2) 
[jac,j1] = jacmat( eq , [x(1),x(2)], x );


for i=1:2
eval(ch(i,:));
end
for i=1:2
eval(eq(i,:));
end
