function [us,tus] = wforu (w,isF,ts,tx,tf)

dim = [];
m=1;
for k=1:isF
    for l=2:length(tx(:,k))
        vector = tx(l,k)-tx(l-1,k);
        if vector == 0
        dim(m)=l-2;
        m=m+1;
        break;
        end
    end
end
n =1;
for i=1:isF
    for j = 1:dim(i)-1;
        for p=n:n+9
        us(p+1,i)=w(j,i);
        tus(p+1,i) = p*tf/((dim(i)-1)*9);
        end
        n = p+1;
    end
    n =1;
end
    us(1,:)=us(2,:);
    tus(1,:) = 0.0;