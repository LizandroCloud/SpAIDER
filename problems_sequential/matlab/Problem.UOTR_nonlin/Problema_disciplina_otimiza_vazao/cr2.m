
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

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

timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
c = [];


timeL=0;
%equality  derivative
dceq = [];


timeL=0;
%inequality  derivative
dc = [];


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
