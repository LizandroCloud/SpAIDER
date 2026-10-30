function [PD,ires]=jac(t,x,xp,cj,rpar)
 global est;
 rpar = rpar';
 ws(:,1)=rpar(3:3+rpar(1,1)-1,1); 
 ks =  rpar(3+rpar(1,1),1); 
 ts(:,1)=rpar(3+rpar(1,1)+1:3+rpar(1,1)+rpar(2,1),1); 
 reg =  rpar(3+rpar(1,1)+1+rpar(2,1),1);
 n_c = rpar(3+rpar(1,1)+rpar(2,1)+2,1);
 n_s = rpar(3+rpar(1,1)+rpar(2,1)+3,1);
 ts_frozen = rpar(3+rpar(1,1)+rpar(2,1) + 4:3+rpar(1,1)+rpar(2,1) + 3+rpar(2,1),1);
 points = rpar(2+rpar(1,1) + rpar(2,1) + 5 + rpar(2,1) :3+rpar(1,1) + rpar(2,1) + 4 + rpar(2,1)+ rpar(length(rpar),1),1);
 kd=3+rpar(1,1) + rpar(2,1) + 4 + rpar(2,1)+ rpar(length(rpar));
 ns_a = rpar(kd:length(rpar)-1,1);


ts_frozen = ts_frozen';
points = points';
k1 = 0.053;  % l/(mol min)
k2 = 0.128; % l/(mol min)
cb_in = 5;   % (mol/l)


% Regularization condition...
for i = 1:length(ws(:,1)) 
if reg==1
     u(i) = u_reg1(t,te(i,:),ws(i,:),ks,ne(i));
else
    if est==0
        u(i) = u_reg0(t,te(i,:),ws(i,:),ks,ne(i));
       else
        u(i) = ws(i);
    end
 end
end


jc1 = [         -k1*x(2)-u/x(3)-cj,              -k1*x(1),           u/x(3)^2*x(1)
              -k1*x(2), -k1*x(1)-4*k2*x(2)-u/x(3)-cj,  -u/x(3)^2*(cb_in-x(2))
                   0,                   0,                   -cj];

PD = jc1;


PD=PD'; % The iteration matrix must be transposed (due to Matlab storage style)
ires = 0;
end