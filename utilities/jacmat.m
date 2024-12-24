function [jac,j1] = jacmat(dfdu , v , x )

for i=1:length(dfdu)
    df(i,:)=jacobian(dfdu(i),v);
end
j1 = df;
for j=1:length(dfdu)
    for i=1:length(v)
        dfduu(j,:)=subs(df( j,: ),v(i),x(i),0);
        df(j,:) = dfduu(j,:);
    end
% jac(j) = sum(df(j,:))';
end
jac = df;
jac=roundn(jac,-2);

end