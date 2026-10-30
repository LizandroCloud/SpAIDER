function [sys,x0,str,ts,simStateCompliance]=transmitter(t,x,u,flag,cinit,par)
%transmitter - transmitter model S-Function
%
%   This model describes the dynamic behavior of a transmitter.
%
%   [sys,x0,str,ts,simStateCompliance]=transmitter(t,x,u,flag,cinit,par)
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
%   par{1}: RangeMin    Transmitter minimum range           []
%   par{2}: RangeMax    Transmitter maximum range           []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): PVinput      Input signal from process variable   []
%   u(2): SPinput      Set-point from user/master loop   	[]
%
%   Outputs:
%   y(1): PVoutput  Bounded process variable output signal  []
%   y(2): SPoutput  Bounded set-point output signal  []
%
%   Initial conditions set:
%   There are no states for this model.

%   2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%----- Model parameters:
RangeMin	= par{1};   %Transmitter minimum range  []
RangeMax	= par{2};	%Transmitter maximum range  []

%----- Inputs:
PVinput   		= u(1);
SPinput  		= u(2);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;   %Number of continuous states
sys(2)=0;   %Number of discrete states
sys(3)=2;   %Number of outputs
sys(4)=2;   %Number of inputs
sys(5)=0;   %Reserved for root finding. Must be zero.
sys(6)=1;   %Direct Feedthrough S-function parameter
sys(7)=1;   %Number of sample times (only for discrete-time stystem)

%----- States initial condition:
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
%----- Model internal variables:
%Bounded process variable output signal to controller input
PVoutput = (PVinput - RangeMin)/(RangeMax - RangeMin);

%Bounded set-point output signal to controller input
SPoutput = (SPinput - RangeMin)/(RangeMax - RangeMin);

%----- Outputs for sys output flag=3:
sys(1)=PVoutput;
sys(2)=SPoutput;
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