function [sys,x0,str,ts,simStateCompliance]=pump_water(t,x,u,flag,cinit,par)
%pump_water - Water brine flow pump non-linear model S-Function
%
%   This model describes the dynamic behavior of a pump with a liquid
%   water brine inlet stream for oil & gas upstream processes, using
%   the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=pump_water(t,x,u,flag,cinit,par)
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
%   par{1}: Cws     Water salt concentration                [g/L]
%   par{2}: h0      Shut-off head                         	[m]
%   par{3}: A       Pressure-flow pump relation coeff.	 	[s^2/m^5]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): Qw_in     Inlet liq. std. vol. flow rate          [Sm^3/d]
%   u(2): P_in      Inlet pressure                          [kPa]
%   u(3): T_in      Inlet temperature                       [K]
%
%   Outputs:
%   y(1): P_out     Outlet pressure                         [kPa]
%
%   Fixed set for steady-state condition:
% 	The model has no dynamics, for which the call with flag=3
% 	computes a steady-state condition.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
Cws = par{1};   %Water salt concentration [g/L]
h0  = par{2};   %Shut-off head [m]
A   = par{3};   %Pressure-flow pump relation coeff. [s^2/m^5]
   
%Water salt mass fraction []
ws=Cws/(999.054+Cws);

%Water brine mass density [kg/m^3]
rWstd=WaterStdMassDens(PP,ws);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=1;%Number of outputs
sys(4)=3;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States at steady-state condition:
x0=[];

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
sys=[];
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,      
%----- Inputs:
Qw_in	= u(1);
P_in    = u(2);
T_in	= u(3);

%----- Model internal variables:
%Inlet water formation volume factor [m^3/Sm^3]
Bw_in=WaterFVF(PP,P_in,T_in);

%Inlet water brine mass density [kg/m^3]
rW_in=rWstd/Bw_in;

%----- Output equations:
%Pressure rise [kPa]
DP=(1e-3)*rW_in*9.81*(h0-A*(Qw_in/(24*3600)*Bw_in)^2);

%Outlet pressure [kPa]
P_out=P_in+DP;

%----- Outputs for sys output flag=3:
sys(1)=P_out;
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