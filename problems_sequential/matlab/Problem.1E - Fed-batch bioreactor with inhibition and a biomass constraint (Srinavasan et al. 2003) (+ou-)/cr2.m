
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

mim = 0.53;
km = 1.2  ;
ki = 22   ;
Yx = 0.4  ;
Yp = 1    ;
v = 0.5   ;
Sin = 20  ;
xo = 1    ;
So = 0    ;
Po = 0    ;
Vo = 2;   ;
xmax= 3   ;
smax= 5   ;

timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,1) - xmax;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,2) - smax;
end



timeL=0;
%equality  derivative
dceq = [];


timeL=0;
%inequality  derivative
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (1)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(1)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:);
end

for i=1:length(df(1,1,:))
end


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
