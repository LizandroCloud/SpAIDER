function [J] = obf1(f,t_ult,t_out,y,stages)
Ad = 2.1428*10^17   ;
Ed = 14897          ;
Ap = 3.816*10^10    ;
Ep = 3557           ;
At = 4.50*10^12     ;
Et = 843            ;
fi = 0.6            ;
DeltaHr = -6.99*10^4;
hA = 1.05*10^6      ;
rhoCp = 1506        ;
rhocCpc = 4043      ;
Mm = 104.14         ;
Qim = 108           ;
Qi = 108            ;
Qmm = 378           ;
Qc =1               ;
V = 3000            ;
Vc = 3312.4         ;
If = 0.5888         ;
Mf = 8.6981         ;
Tf = 330            ;
Tcf = 295           ;

J=- 1*f(end,8);
end