function [sys,x0,str,ts,simStateCompliance]=valv_water(t,x,u,flag,cinit,par)
%valv_water - Watre brine flow valve non-linear model S-Function
%
%   This model describes the dynamic behavior of a valve with a water
%   brine only inlet stream for oil & gas upstream processes, using
%   the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=valv_water(t,x,u,flag,cinit,par)
%
%   sys is the S-function output, which depends on which flag is 
%   being inputted by simulink. x0 are the states initial condi-
%   tions. See the S-function Level-1 MATLAB reference page for more
%   details about str, ts and simStateCompliance S-Function ouputs.
%
%   t is the independent variable time. x are the state variables.
%   u are the input variables, combined into one vector of manipula-
%   ted and disturbs. flag is the Simulink variable inputted to get a
%   specific output from the S-function. cinit are variable set for 
%	state initial condition computation. par are the model parameters.
%
%   Model parameters:
%   par{1}: Cws     Water salt concentration                	[g/L]
%   par{2}: Cv      Valve coefficient constant US units      	[gpm]
%   par{3}: FL      Liquid pressure recovery factor          	[]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): P_in      Inlet pressure                           [kPa]
%   u(2): T_in      Inlet temperature                        [K]
%   u(3): P_out     Outlet pressure                          [kPa]
%   u(4): c     	Valve stem position            			 []
%
%   Outputs:
%   y(1): Qw        Water std. vol. flow rate                [Sm^3/d]
%   y(2): T         Temperature                              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
Cws         = par{1};   %Water salt concentration [g/L]
Cv          = par{2};   %Valve constant US units [gpm]
FL          = par{3};   %Liquid pressure recovery factor []

%Water salt mass fraction []
ws=Cws/(999.054+Cws);

%Water brine mass density [kg/m^3]
rWstd=WaterStdMassDens(PP,ws);

%Valve coefficient constant [m^2]
Av=Cv*(0.0000630902)*sqrt(PP.rho_H2Ostd*14.696/101325);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=4;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States initial condition:
x0=[];

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
%----- State derivatives for sys output flag=1
sys=[];
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Inputs:
P_in    = u(1);
T_in	= u(2);
P_out   = u(3);
c       = u(4);

%----- Model internal variables:
%Pressure drop [Kpa]
DP=P_in-P_out;

%Liquid phase vapour pressure [kPa]
PvW=(1e2)*10^(5.11564-1687.537/((T_in-273.15)+230.17));

%Molar mixing rule for liquid phase pseudo-critical pressure [kPa]
PpcW=22064;

%Liquid critical pressure ratio factor
FF=0.5*(1+tanh((PvW-PpcW)/(100)/1e-6))*(0.96-0.28) +...
 0.5*(1+tanh((PpcW-PvW)/(100)/1e-6))*(0.96-0.28*sqrt(PvW/PpcW));

%Linear intrinsic characteristic
fc=c;

%Pressure drop at choked flow
DPchok=FL^2*(P_in-FF*PvW);

%----- Output equations:
%Pressure-flow relation
DPflow=0.5*(1+tanh(((DPchok-DP)+1e-3)/1e-4))*DP...
 +0.5*(1+tanh((DP-DPchok)/1e-4))*DPchok;
Qw=(24*3600)/rWstd*fc*Av*sqrt(rWstd*DPflow*1e3);

%Energy balance: isothermal flow
T=T_in;

%----- Outputs for sys output flag=3:
sys(1)=Qw;
sys(2)=T;
sys=sys(:);
%--------------------------------------------------------------------
%----- S-FUNCTION OTHER FLAGS ---------------------------------------
    case 4,
sys=[];
    case 9,
sys=[];
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end