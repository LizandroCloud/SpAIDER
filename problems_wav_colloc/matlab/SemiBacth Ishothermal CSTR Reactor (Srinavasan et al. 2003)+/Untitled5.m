J=2;
t=1;

for j = 0:J
    for k = 0:(2^j - 1);
        n = k + 2^j + 1;
        hi=haar2(n,j,k,t);
        h(n)=hi;
        hi=0.0;
    end
end