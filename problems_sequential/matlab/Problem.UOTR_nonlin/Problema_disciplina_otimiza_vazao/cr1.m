function [ceq c] = cr1(f,t_ult,t_out,u,y,stages)
% r1amax=0.85              ;
% r1bmax=0.022             ;
% r2mmax=0.018             ;
% r3max=0.09               ;
% r4max=0.03               ;
% r5max=0.03               ;
% r6max=0.75               ;
% r7max=0.03               ;
% r7mmax=0.012             ;
% k1a=0.08                 ;
% k1b=0.001                ;
% k2m=0.009                ;
% k3=0.09                  ;
% k4=0.05                  ;
% k5=0.05                  ;
% k6=0.05                  ;
% k7=0.001                 ;
% k7m=0.03                 ;
% y1=1                     ;
% y2m=1                    ;
% y3=0.68                  ;
% y4=0.7                   ;
% y5=1                     ;
% y6=0.7                   ;
% y7=1.1                   ;
% ke1a=20*k1a              ;
% remax=0.001              ;
% beta=0.1                 ;
% rgmaxchapeu=0.0575       ;
% kglinha=0.01             ;
% kg=0.05                  ;
% D = 0.031                ;
% ycm1=0.2                 ;
% ycm2=0.15                ;
% ycm3=0.15                ;
% ycm4=(1-ycm1-ycm2-ycm3)/2;
% s1f=0.95                 ;
% s2f=0.4238               ;
% 
% timeL=0;
% 
% %equality
% ceq = [];
% 
% timeL=0;
% 
% %inequality
% c(1) = 0.0 - f(end,1);
% c(2) = f(end,1)-0.3;
% c(3) = 0.0 - f(end,2);
% c(4) = f(end,2) -0.25;
% c(5) = 0.0 - f(end,3);
% c(6) = f(end,3)-0.1;
% c(7) = 0.05 - f(end,4);
% c(8) = f(end,4) - 0.2;
% c(9) = 0.3 - f(end,5);
% c(10) = f(end,5) - 0.6;
% c(11) = 0.0 - f(end,6);
% c(12) = f(end,6) - 0.25;
% c(13) = 0.24 - f(end,7);
% c(14) = f(end,7) - 0.5;
% c(15) = 0.20 - f(end,8);
% c(16) = f(end,8)- 0.4;
% c(17) = 0.40 - f(end,9);
% c(18) = f(end,9)- 0.75;
% c(19) = 0.0- f(end,10);
% c(20) = f(end,10)- 0.005;
% c(21) = 0.0 - f(end,11);
% c(22) = f(end,11)- 0.0;
% c(23) = 0.005-f(end,12);
% c(24) = f(end,12)- 0.01;
% c(25) = 0.0-f(end,13);
% c(26) = f(end,13)- 0.0;
% c(27) = 0.002-f(end,14);
% c(28) = f(end,14)- 0.005;
% c(29) = 0.0-f(end,15);
% c(30) = f(end,15)- 0.0;
% for time=1:40
%     c(29+time) = f(time,7) - 0.25;
% end
% ki=time;
% for time=1:40
%     c(29+time) = f(time,6) - 0.2;
% end
c=[];
ceq=[];



% c(3) = f(:,9)-1.25*ones(40,1);

end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
