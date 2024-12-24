function pzo = f_zo(ks,z0,x0_sm)

if ks==1
    x = cell2mat(z0);
    pzo = x(ks,:);
else
    x = cell2mat(x0_sm);
    pzo = x(ks-1,:);
end

end