
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)


timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,1)-40;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,3)-25;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,4)-10;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,1)+0 ;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,3)+0 ;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,4)+0 ;
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
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(1)*df(time,3,:)+(0)*df(time,4,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(1)*df(time,4,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (1)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(0)*df(time,4,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(1)*df(time,3,:)+(0)*df(time,4,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(0)*df(time,3,:)+(1)*df(time,4,:);
end

for i=1:length(df(1,1,:))
end


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
