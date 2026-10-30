function sys=sfunc2ss(sfunc,cinit,u0,par)
%sfunc2ss - S-Function nonlinear model to state-space linear time in-
%variant (LTI) model conversion.
%
%   sys = sfunc2ss(sfunc,x0,u0,cinit,par)
%
%   sys is the state-space linear model from sys=ss(A,B,C,D). A,B,C,D
%   are the state-space matrices obtained from the jacobian matrices
%   of Level-1 MATLAB S-function sfunc.
%
%   sfunc is the S-function function handle. cinit,u0 are the states and
%	inputs condition of the linearization point. par are the S-function 
%	model parameters.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

[sys,x0]=sfunc([],[],[],0,cinit,par);
Nx=sys(1);
Nu=sys(4);

if Nx>0
    fhand=@(t,x)sfunc(t,x,u0,1,cinit,par);
    fval=fhand(0,x0);
    [A,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));
else
    A=[];
end

if Nx>0 && Nu>0
    fhand=@(t,u)sfunc(t,x0,u,1,cinit,par);
    fval=fhand(0,u0);
    [B,~]=numjac(fhand,0,u0,fval,zeros(Nu,1));
else
    B=[];
end

if Nx>0
    fhand=@(t,x)sfunc(t,x,u0,3,cinit,par);
    fval=fhand(0,x0);
    [C,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));
else
    C=[];
end

if Nu>0
    fhand=@(t,u)sfunc(t,x0,u,3,cinit,par);
    fval=fhand(0,u0);    
    [D,~]=numjac(fhand,0,u0,fval,zeros(Nu,1));
else
    D=[];
end

sys=ss(A,B,C,D);
end