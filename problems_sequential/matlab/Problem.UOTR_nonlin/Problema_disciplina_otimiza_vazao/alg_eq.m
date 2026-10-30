function [yout] = alg_eq(x,u)
r1amax=0.85              ;
r1bmax=0.022             ;
r2mmax=0.018             ;
r3max=0.09               ;
r4max=0.03               ;
r5max=0.03               ;
r6max=0.75               ;
r7max=0.03               ;
r7mmax=0.012             ;
k1a=0.08                 ;
k1b=0.001                ;
k2m=0.009                ;
k3=0.09                  ;
k4=0.05                  ;
k5=0.05                  ;
k6=0.05                  ;
k7=0.001                 ;
k7m=0.03                 ;
y1=1                     ;
y2m=1                    ;
y3=0.68                  ;
y4=0.7                   ;
y5=1                     ;
y6=0.7                   ;
y7=1.1                   ;
ke1a=20*k1a              ;
remax=0.001              ;
beta=0.1                 ;
rgmaxchapeu=0.0575       ;
kglinha=0.01             ;
kg=0.05                  ;
D = 0.031                ;
ycm1=0.2                 ;
ycm2=0.15                ;
ycm3=0.15                ;
ycm4=(1-ycm1-ycm2-ycm3)/2;
s1f=0.95                 ;
s2f=0.4238               ;


y(1)=(r1amax*(x(1)/(k1a+x(1)+((x(1)^2)/k1a)))*(x(7)/(x(7)+kglinha))*(x(10)/(((10^-7)+remax)/(rgmaxchapeu+beta)))                 );

y(2)=(r1bmax*(x(1)/(k1b+x(1)+((x(1)^2)/k1b)))*(x(7)/(x(7)+kglinha))*(x(11)/(((10^-7)+remax)/(rgmaxchapeu+beta)))                 );

y(3)=(r3max*(x(2)/(k3+x(2)+((x(2)^2)/k3)))*(x(7)/(x(7)+kglinha))*(x(12)/(((10^-7)+remax)/(rgmaxchapeu+beta)))                    );

y(4)=(r7max*(x(2)/(k7+x(2)+((x(2)^2)/k7)))*(x(7)/(x(7)+kglinha))*(x(15)/(((10^-7)+remax)/(rgmaxchapeu+beta)))                    );

y(5)=(u(1)                                                                                                                       );

y(6)=(r6max*(x(3)/(k6+x(3)))*(x(7)/(x(7)+kglinha))*(x(14)/(((10^-7)+remax)/(rgmaxchapeu+beta)))                                  );

y(7)=(r5max*(x(4)/(k5+x(4)))*(x(7)/(x(7)+kglinha))                                                                               );

y(8)=(r2mmax*(x(1)/(k2m+x(1)))*(x(7)/(x(7)+kglinha))                                                                             );

y(9)=(r7mmax*(x(2)/(k7m+x(2)))*(x(7)/(x(7)+kglinha))                                                                             );

y(10)=(remax*(x(1)/(ke1a+x(1)))                                                                                                   );

y(11)=(remax*(x(1)/(k1b+x(1)))                                                                                                    );

y(12)=(remax*(x(2)/(k3+x(2)))                                                                                                     );

y(13)=(remax*(x(3)/(k4+x(3)))                                                                                                     );

y(14)=(remax*(x(3)/(k6+x(3)))                                                                                                     );

y(15)=(remax*(x(2)/(k7+x(2)))                                                                                                     );

y(16)=(rgmaxchapeu*(x(3)/(x(3)+kg))*(x(4)/(x(4)+kg))*(x(5)/(x(5)+kg))*(x(6)/(x(6)+kg))                                            );

y(17)=(y(1)/y(2)                                                                                                                  );

y(18)=(y(2)/y(2)                                                                                                                  );

y(19)=(y(3)/ y(4)                                                                                                                 );

y(20)=(y(4)/y(3)                                                                                                                  );

y(21)=(y3*y(3)/y3*y(3)                                                                                                            );

y(22)=(y4*y(5)/y4*y(5)                                                                                                            );

y(23)=(y1*(y(1)*y(17)+y(2)*y(18))+y2m*y(8)+y3*y(3)*y(19)+y7*(y(4)*y(20)+y(9))+(y4-1)*y(7)-y(6)*y(17)+y(16)*(1-ycm1-ycm2-ycm3-ycm4));

y(24)=(y(1)/(y(1)+y(2))                                                                                                           );

y(25)=(y(2)/(y(1)+y(2))                                                                                                           );

y(26)=(y(3)/(y(3)+y(4))                                                                                                           );

y(27)=(y(4)/(y(3)+y(4))                                                                                                           );

y(28)=(y3*y(3)/(y3*y(3)+y4*y(5))                                                                                                  );

y(29)=(y4*y(5)/(y3*y(3)+y4*y(5))                                                                                                  );

yout=y;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
