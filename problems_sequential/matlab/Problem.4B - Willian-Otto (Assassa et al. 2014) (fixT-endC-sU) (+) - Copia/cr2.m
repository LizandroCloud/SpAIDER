
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

UA = 1560     ;
V = 9450      ;
Vj = 2000     ;
dh = -69919.56;
Rg = 1.987    ;
roj = 1       ;
ros = 0.915   ;
cpj = 4045.7  ;
cps = 1647.27 ;
T0 = 333      ;
Tj0 = 294     ;
ci0 = 0.981   ;
cm0 = 8.63    ;
cb0 = 1.05    ;

timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,8) - 380 ;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = -f(time,8) + 320;
end



timeL=0;
%equality  derivative
dceq = [];


timeL=0;
%inequality  derivative
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:)+(0)*df(time,5,:)+(0)*df(time,6,:)+(0)*df(time,7,:)+(1)*df(time,8,:)+(0)*df(time,9,:)+(0)*df(time,10,:) ;
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:)+(0)*df(time,5,:)+(0)*df(time,6,:)+(0)*df(time,7,:)+(-1)*df(time,8,:)+(0)*df(time,9,:)+(0)*df(time,10,:);
end

for i=1:length(df(1,1,:))
end


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
