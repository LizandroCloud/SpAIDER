%Script to run PID_incr Level-1 MATLAB S-Function
clear

%----- Setting model parameters:
action = {'Reverse'};

%----- Linearization point:
mode		= 0;	%Controller mode: Manual (1), Auto (0)  []
PV  		= 0.9;	%Bounded input signal from transmitter  []
SP  		= 1.0;	%Bounded set-point signal from trans.	[]
KC  		= 2.0;	%Proportional gain                     	[]
TI  		= 30;	%Integral time        					[s]
TD  		= 2.0;	%Derivative time        				[s]

unboutp 	= 0.5;		%Unbounded output  						[]
errorTrack 	= SP-PV;	%Error tracking       					[]
dderivFiltdt= 0;		%Deriv. action w/filter for incr. PID	[]
MVinput 	= unboutp;	%MV tracking       						[]

%----- dynamic simulation
x0 = [...
	unboutp; ...
	errorTrack; ...
	dderivFiltdt; ...
	MVinput];

u = [...
	mode; ...
	PV; ...  
	SP; ...  
	KC; ...  
	TI; ...  
	TD];

tspan=0:1:300;
[t,y]=ode15s(@(t,x)PID_incr(t,x,u,1,[],action),tspan,x0);
plot(t,y);