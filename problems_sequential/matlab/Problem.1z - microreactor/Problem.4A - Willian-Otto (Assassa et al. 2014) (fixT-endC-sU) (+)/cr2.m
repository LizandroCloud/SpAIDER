
function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)

E1 = -6666.7  ;
E2 = -8333.3  ;
E3 = -11111   ;
k1 = 1.6599e6 ;
k2 = 7.2117e8 ;
k3 = 2.6745e12;
dh1 = -263.8  ;
dh2 = -158.3  ;
dh3 = -226.3  ;
cp = 4.184    ;
cpx = 5554.1  ;
cex = 125.91  ;
Tin = 35      ;
Aref = 9.2903 ;
mref = 2.1052 ;
U = 2.3082e-1 ;
Ma = 100      ;
Mb = 100      ;
Mc = 200      ;
Me = 200      ;
Mg = 300      ;
Mp = 100      ;
roa = 1       ;
rob = 0       ;
roc = 0       ;
roe = 0       ;
rog = 0       ;
rop = 0       ;
tro = 65      ;
Vo = 2        ;
Fscal = 1000  ;

timeL=0;

%equality
ceq = [];

timeL=0;

%inequality
c(1) = f(end,8) - 5.0;


timeL=0;
%equality  derivative
dceq = [];


timeL=0;
%inequality  derivative
dc(1,:) = (0)*df(end,1,:)+(0)*df(end,2,:)+(0)*df(end,3,:)+(0)*df(end,4,:)+(0)*df(end,5,:)+(0)*df(end,6,:)+(0)*df(end,7,:)+(1.0)*df(end,8,:)+(0)*df(end,9,:);


dc=transp(dc);
dceq=transp(dceq);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
