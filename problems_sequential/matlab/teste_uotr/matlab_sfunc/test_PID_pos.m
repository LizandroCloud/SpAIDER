%Script to run uvrc Level-1 MATLAB S-Function
clear

%----- Setting model parameters:
action = {'Reverse'};

%----- Linearization point:
mode		= 0;	%Controller mode: Manual (1), Auto (0)  []
PV  		= 0.9;	%Bounded input signal from transmitter  []
SP  		= 1.0;	%Bounded set-point signal from trans.	[]
KC  		= 1.0;	%Proportional gain                     	[]
TI  		= 60;	%Integral time        					[s]
TD  		= 2.0;	%Derivative time        				[s]

AWBTintterm = 0.5;			%Integral action with AWBT  	[]
derivFilt 	= SP-PV;		%Derivative action with filter  []
MVinput 	= AWBTintterm;	%MV tracking       				[]

%----- dynamic simulation
x0 = [...
	AWBTintterm; ...
	derivFilt; ...
	MVinput];

u = [...
	mode; ...
	PV; ...  
	SP; ...  
	KC; ...  
	TI; ...  
	TD];

tspan=0:1:300;
[t,y]=ode15s(@(t,x)PID_pos(t,x,u,1,[],action),tspan,x0);
plot(t,y);