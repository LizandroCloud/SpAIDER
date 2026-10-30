function [ceq c] = cr1(f,t_ult,t_out,u,y,stages, ns_a)

% constraints without derivatives; same values as cr2
cb_in = 5      ;
Vo = 1         ;
cao=0.72       ;
cbo=0.05       ;
cd_mais = 0.15 ;
cb_mais = 0.025;

%equality
ceq = [];

%inequality
c(1) = ((0.5/f(end,3))*((f(end,1)+cb_in-f(end,2))*f(end,3) - (cao+cb_in-cbo)*Vo)) - cd_mais;
c(2) = f(end,2) - cb_mais;
end
