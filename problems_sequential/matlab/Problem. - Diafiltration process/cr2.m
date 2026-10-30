
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

s1 = 68.1250*10^-9 ;
s2 = -56.4512*10^-6;
s3 = 32.5553*10^-3 ;
s4 = -4.3529*10^-9 ;
s5 = 3.3216*10^-6  ;
s6 = -2.7141*10^-3 ;
z1 = 7.8407*10^-6  ;
z2 =  -4.0507*10^-3;
z3 = 1.0585        ;
z4 =  1.2318*10^-9 ;
z5 = -9.7660*10^-6 ;
z6 = -1.1677*10^-3 ;
w1 = -0.0769*10^-6 ;
w2 =  -0.0035*10^-3;
w3 = 0.0349*10^-3  ;
w4 =  0.9961       ;

timeL=0;

%equality
ceq(1) =  f(end,3) - 0.01 ;
timeL=0;

%inequality
for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = -f(time,3) + 0.01 ;
end

timeL=length(f(:,1));;


for time=1+0+0:length(f(:,1))+0+0
c(time+timeL) = f(time,3) - 0.035 ;
end



timeL=0;
%equality  derivative
dceq(1,:) = (0)*df(end,1,:)+(0)*df(end,2,:)+(1.0)*df(end,3,:);


timeL=0;
%inequality  derivative
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(-1.0)*df(time,3,:);
end
timeL=length(f(:,1));
for time=1+0+0:length(f(:,1))+0+0
dc(time+timeL,:) = (0)*df(time,1,:)+(0)*df(time,2,:)+(1.0)*df(time,3,:) ;
end

for i=1:length(df(1,1,:))
end


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
