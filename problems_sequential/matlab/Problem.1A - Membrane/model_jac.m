function [res,ires] = model_jac(t,x,xp,rpar)
global DEBUG;
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

res(1)= ( -k1*x(1)*x(2)-( u/x(3) )*x(1)) - xp(1);
res(2)= ( -k1*x(1)*x(2) - 2*k2*x(2)^2 + ( (u/x(3))*(cb_in - x(2)))) - xp(2); 
res(3)= u - xp(3);


jc1 = [         -k1*x(2)-u/x(3),              -k1*x(1),           u/x(3)^2*x(1),            -1/x(3)*x(1)
              -k1*x(2), -k1*x(1)-4*k2*x(2)-u/x(3),  -u/x(3)^2*(cb_in-x(2)),     1/x(3)*(cb_in-x(2))
                   0,                   0,                   0,                   1];

jk=double(jc1);
j0=0;
jk(:,1:n_s)=jc1(:,1:n_s);
jk(:,n_s+1:n_s+n_c)=0;



for j=1:n_c
    
%      tx = [t ts_frozen(1,2:ns_a(j)+1)];
    tx = [t ts_frozen(points(j,1:length(points)))];
    tx = sort(tx);
    check = find(tx == t,1,'last');
    
    for is = 1:ns_a(j) %length(ws)
        if ((ks)== is)
             uws = 1;
        else   % do contrario...
             uws = 0;
        end
                jk(:,n_s+j)=uws*jc1(:,n_s+j);
                k=1:(n_s);
                res(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0)+ jk(k,n_s+j) - xp(n_s*is+k+j0);
    end 
j0=length(res)-n_s;
end




ires = 0;

if ( DEBUG > 9 )
  buf1 = sprintf('%.4e ',y);
  buf2 = sprintf('%.4e ',yp);
  buf3 = sprintf('%.4e ',res);
  disp(sprintf('(Debug) dydt.m: t=%.2e y=%s yp=%s res=%s', ...
               t,buf1,buf2,buf3));
  if ( DEBUG > 99 )
    cbuf = input(sprintf('\tRETURN to continue : '));
  end

end
