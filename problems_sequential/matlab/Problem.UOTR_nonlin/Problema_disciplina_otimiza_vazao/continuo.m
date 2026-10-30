%Continuous Reactor

function ds = continuo(t,x,par)

r1amax=par(1);
r1bmax=par(2);
r2mmax=par(3);
r3max=par(4);
r4max=par(5);
r5max=par(6);
r6max=par(7);
r7max=par(8);
r7mmax=par(9);

k1a=par(10);
k1b=par(11);
k2m=par(12);
k3=par(13);
k4=par(14);
k5=par(15);
k6=par(16);
k7=par(17);
k7m=par(18);

y1=par(19);
% y1b=par(20);
y2m=par(20);
y3=par(21);
y4=par(22);
y5=par(23);
y6=par(24);
y7=par(25);
% y7m=par(27);

ke1a=par(26);

remax=par(27);
beta=par(28);

rgmaxchapeu=par(29);
kglinha=par(30);
kg=par(31);
D=par(32);

ycm1=par(33);
ycm2=par(34);
ycm3=par(35);
ycm4=par(36);

s1f=par(37);
s2f=par(38);

%x(1)=s1
%x(2)=s2
%x(3)=m1
%x(4)=m2
%x(5)=m3
%x(6)=m4
%x(7)=c'
%x(8)=x
%x(9)=L
%x(10)=e1a
%x(11)=e1b
%x(12)=e3
%x(13)=e4
%x(14)=e6
%x(15)=e7

r1a=r1amax.*(x(1)./(k1a+x(1)+((x(1).^2)./k1a))).*(x(7)./(x(7)+kglinha)).*(x(10)./(((10^-7)+remax)./(rgmaxchapeu+beta)));
r1b=r1bmax.*(x(1)./(k1b+x(1)+((x(1).^2)./k1b))).*(x(7)./(x(7)+kglinha)).*(x(11)./(((10^-7)+remax)./(rgmaxchapeu+beta)));
r3=r3max.*(x(2)./(k3+x(2)+((x(2).^2)./k3))).*(x(7)./(x(7)+kglinha)).*(x(12)./(((10^-7)+remax)./(rgmaxchapeu+beta)));
r7=r7max.*(x(2)./(k7+x(2)+((x(2).^2)./k7))).*(x(7)./(x(7)+kglinha)).*(x(15)./(((10^-7)+remax)./(rgmaxchapeu+beta)));

r4=r4max.*(x(3)./(k4+x(3))).*(x(7)./(x(7)+kglinha)).*(x(13)./(((10^-7)+remax)./(rgmaxchapeu+beta)));
r6=r6max.*(x(3)./(k6+x(3))).*(x(7)./(x(7)+kglinha)).*(x(14)./(((10^-7)+remax)./(rgmaxchapeu+beta)));

r5=r5max.*(x(4)./(k5+x(4))).*(x(7)./(x(7)+kglinha));

r2m=r2mmax.*(x(1)./(k2m+x(1))).*(x(7)./(x(7)+kglinha));
r7m=r7mmax.*(x(2)./(k7m+x(2))).*(x(7)./(x(7)+kglinha));

re1a=remax.*(x(1)./(ke1a+x(1)));
re1b=remax.*(x(1)./(k1b+x(1)));
re3=remax.*(x(2)./(k3+x(2)));
re4=remax.*(x(3)./(k4+x(3)));
re6=remax.*(x(3)./(k6+x(3)));
re7=remax.*(x(2)./(k7+x(2)));


rgchapeu=rgmaxchapeu.*(x(3)./(x(3)+kg)).*(x(4)./(x(4)+kg)).*(x(5)./(x(5)+kg)).*(x(6)./(x(6)+kg));

vc1a=r1a./(max(r1a,r1b));
vc1b=r1b./(max(r1a,r1b));
vd3=r3./(max(r3,r7));
vd7=r7./(max(r3,r7));
vc3=y3.*r3./(max(y3*r3,y4*r4));
vc4=y4.*r4./(max(y3*r3,y4*r4));


rg=y1.*(r1a.*vc1a+r1b.*vc1b)+y2m.*r2m+y3.*r3.*vd3+y7.*(r7.*vd7+r7m)+(y4-1).*r5-r6.*vc1a+rgchapeu.*(1-ycm1-ycm2-ycm3-ycm4);

uc1a=r1a./(r1a+r1b);
uc1b=r1b./(r1a+r1b);
ud3=r3./(r3+r7);
ud7=r7./(r3+r7);
uc3=y3.*r3./(y3.*r3+y4.*r4);
uc4=y4.*r4./(y3.*r3+y4.*r4);


ds(1)=-(r1a.*vc1a + r1b.*vc1b +r2m).*x(8) +D.*(s1f-x(1)); %ds1/dt
ds(2)=-(r3.*vc3.*vd3+r7.*vd7+r7m).*x(8)+D.*(s2f-x(2)); %ds2/dt
ds(3)=y1.*(r1a.*vc1a+r1b.*vc1b)+y5.*r5-r4.*vc4-r6.*vc1a-ycm1.*rgchapeu-rg.*x(3); %dm1/dt
ds(4)=y3.*r3.*vc3.*vd3+y4.*r4.*vc4-r5-ycm2.*rgchapeu-rg.*x(4); %dm2/dt
ds(5)=y2m.*r2m-ycm3.*rgchapeu-rg.*x(5); %dm3/dt
ds(6)=y7.*(r7.*vd7+r7m)-ycm4.*rgchapeu-rg.*x(6); %dm4/dt
ds(7)=rgchapeu-rg.*x(7); %dc'/dt
ds(8)=(rg-D).*x(8); %dx/dt
ds(9)=y6.*r6.*vc1a.*x(8)-D.*x(9); %dL/dt
ds(10)=(10^-7)+re1a.*uc1a-(rg+beta).*x(10); %de1a/dt
ds(11)=(10^-7)+re1b.*uc1b-(rg+beta).*x(11); %de1b/dt
ds(12)=(10^-7)+re3.*uc3.*ud3-(rg+beta).*x(12); %de3/dt
ds(13)=(10^-7)+re4.*uc4-(rg+beta).*x(13); %de4/dt
ds(14)=(10^-7)+re6.*uc1a-(rg+beta).*x(14); %de6/dt
ds(15)=(10^-7)+re7.*ud7-(rg+beta).*x(15); %de7/dt


ds=ds';




