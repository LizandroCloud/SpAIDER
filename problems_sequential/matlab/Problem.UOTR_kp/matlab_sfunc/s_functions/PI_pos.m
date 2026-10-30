function [sys,x0,str,ts,simStateCompliance]=PI_pos(t,x,u,flag,cinit,par)
%PID_pos - Positional PI controller model S-Function
%
%   This model describes the dynamic behavior of a positional PID
%	controller.
%
%   [sys,x0,str,ts,simStateCompliance]=PID_pos(t,x,u,flag,cinit,par)
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
%   par{1}: action      Controller action                	[]
%
%   States:
%   x(1): AWBTintterm   Integral action with AWBT  			[]
%
%   Inputs:
%   u(1): mode      Controller mode: Manual (1), Auto (0)   []
%   u(2): PV      	Bounded input signal from transmitter   []
%   u(3): SP      	Bounded set-point signal from trans.	[]
%   u(4): KC    	Proportional gain                     	[]
%   u(5): TI     	Integral time        					[s]
%   u(6): TD     	Derivative time        					[s]
%
%   Outputs:
%   y(1): MVoutput  Bounded manipulated var. output signal  []
%
%   Initial conditions set:
%   %cinit = x

%   2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%----- Model parameters:
h           = 1e-2;     %regularization function parameter []
action      = par{1};   %Controller action []

%Controller direction
switch action,
	case 'Direct',
	PID_action = -1.0;
		
	case 'Reverse',
	PID_action = 1.0;
end

%----- States:	
AWBTintterm	= x(1);
%derivFilt   = x(2);
%MVinput   	= x(3);

%----- Inputs:
mode   		= u(1);
PV  		= u(2);
SP  		= u(3);
KC			= u(4);
TI			= u(5);
TD			= u(6);

%----- Model internal variables:
%Controller error
error = SP - PV;

%Conventional x decoupled PID law: integral action parameter
KI = KC/TI;

%Conventional x decoupled PID law: derivative action parameter
KD = KC*TD;

%Derivative term filter
%dderivFiltdt = (0.1*(TD+1e-3))\(error - derivFilt);
%dderivFiltdt = (0.1*(TD+1e-2))\(error - derivFilt);

%PID positional law
%unboutp = PID_action*(KC*error + AWBTintterm + KD*dderivFiltdt);
unboutp = PID_action*(KC*error + AWBTintterm);

%Bounded output
%{
MVoutput = ...
	0.5*(1+tanh((2*mode-1)/1e-6))*MVinput ...
	+ 0.5*(1+tanh((1-2*mode)/1e-6))*(...
		0.5*(1+tanh((unboutp-1e-5)/1e-4))*0.5*(1+tanh((0.99999-unboutp)/1e-4))*unboutp ...
		+ 0.5*(1+tanh((unboutp-0.99999)/1e-4)) ...
	);
%}
%{
MVoutput = (...
		0.5*(1+tanh((unboutp-1e-5)/1e-4))*0.5*(1+tanh((0.99999-unboutp)/1e-4))*unboutp ...
		+ 0.5*(1+tanh((unboutp-0.99999)/1e-4)) ...
	);
%}
MVoutput = (...
		0.5*(1+tanh((unboutp-1e-5)/h))*0.5*(1+tanh((0.99999-unboutp)/h))*unboutp ...
		+ 0.5*(1+tanh((unboutp-0.99999)/h)) ...
	);
%----- State equations:
%Integral term with anti-windup and bumpless transfer
%{
dAWBTinttermdt = KI*error* ...
	0.5*(1+tanh((unboutp-1e-5)/1e-6))*0.5*(1+tanh((0.99999-unboutp)/1e-6));
%}

dAWBTinttermdt = KI*error* ...
	0.5*(1+tanh((unboutp-1e-5)/h))*0.5*(1+tanh((0.99999-unboutp)/h));

%MV tracking
%dMVinputdt = 1e-3\(MVoutput - MVinput);
%dMVinputdt = 1e-2\(MVoutput - MVinput);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=1;   %Number of continuous states
sys(2)=0;   %Number of discrete states
sys(3)=1;   %Number of outputs
sys(4)=6;   %Number of inputs
sys(5)=0;   %Reserved for root finding. Must be zero.
sys(6)=1;   %Direct Feedthrough S-function parameter
sys(7)=1;   %Number of sample times (only for discrete-time stystem)

%----- States initial condition:
x0=cinit;
x0=x0(:);

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
%----- State derivatives for sys output flag=1
sys(1)=dAWBTinttermdt;
%sys(2)=dderivFiltdt;
%sys(3)=dMVinputdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Outputs for sys output flag=3:
sys(1)=MVoutput;
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