
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

Ad = 2.1428*10^17   ;
Ed = 14897          ;
Ap = 3.816*10^10    ;
Ep = 3557           ;
At = 4.50*10^12     ;
Et = 843            ;
fi = 0.6            ;
DeltaHr = -6.99*10^4;
hA = 1.05*10^6      ;
rhoCp = 1506        ;
rhocCpc = 4043      ;
Mm = 104.14         ;
Qim = 108           ;
Qi = 108            ;
Qmm = 378           ;
Qc =1               ;
V = 3000            ;
Vc = 3312.4         ;
If = 0.5888         ;
Mf = 8.6981         ;
Tf = 330            ;
Tcf = 295           ;

timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,3)-387                                  ;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = 104.14*(f(time,7)*f(time,5)/(f(time,6)^2))-1.51;
end



timeL=0;
%equality  derivative
dceq = [];


timeL=0;
%inequality  derivative
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(1)*df(time,3,:)+(0)*df(time,4,:)+(0)*df(time,5,:)+(0)*df(time,6,:)+(0)*df(time,7,:)+(0)*df(time,8,:)                                                               ;
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:)+((104.14*f(7))/f(6)^2)*df(time,5,:)+(-(208.28*f(5)*f(7))/f(6)^3)*df(time,6,:)+((104.14*f(5))/f(6)^2)*df(time,7,:)+(0)*df(time,8,:);
end

for i=1:length(df(1,1,:))
end


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
