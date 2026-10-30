function [ceq c] = cr1(f,t_ult,t_out,u,y,stages)

jmax= length(f(1,:));       %states
nk	= length(stages)-1;     %number of elements
ne	= (length(f(:,1)) - 1)   /nk;	%length of the element

%----- Equality constraints ---------------------------------------
%Equality constraint for terminal conditions of states derivatives
ceqdfterm=[...
    f(end,1:end-1)-f(end-1,1:end-1),			... from x(1) to x(27)
	f(end,end)-f(end-1,end) - 1000/(3600*24)];	%dVg,flare/dt = 1000 Mm^3/d
ceqdfterm=ceqdfterm';

%Equality constraint for terminal conditions of controlled outputs
%idy=[1,4,8,11,12,14,16,18,21,23,26,27];
idy=[1,4,8,12,14,16,18,26];
ceqyterm=f(end,idy)-f(1,idy);
ceqyterm=ceqyterm';

%Equality constraint for final and initial conditions between nodal points for multiple-shooting approach
ceq=[ceqdfterm; ceqyterm];
% ko=numel(ceq);
% for j=1:jmax
% 	for i=1:nk-1
% 		ceqNP(i+ko)=f(i*ne+1,j)-f(i*ne,j);
% 	end
% 	ko=ko+i;
% end
% ceqNP=ceqNP(:);
% 
% %outputting equality constraints
% ceq=[ceq; ceqNP];
%--------------------------------------------------------------------
%----- Inequality constraints ---------------------------------------
c=[]; ko=numel(c);

%States upper and lower bounds
fLB = [...
	0.1;...		%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0;...		%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	1;...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1;...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0;...		%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	1600;...	%x(6):  SG_001.P        			Pressure                            [kPa]
	270;...		%x(7):  SG_001.T        			Temperature                         [K]
	1;...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0;...		%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	700;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	340;...		%x(11): V_TO_001.T      			Temperature                         [K]
	1;...		%x(12): TO_001.hW      				Water phase level                   [m]
	340;...		%x(13): TO_001.T       				Temperature                         [K]
	1;...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0;...		%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	150;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	340;...		%x(17): V_TO_002.T      			Temperature                         [K]
	1;...		%x(18): TO_002.hW      				Water phase level                   [m]
	340;...		%x(19): TO_002.T       				Temperature                         [K]
	150;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	290;...		%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	700;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	290;...		%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	0;...		%x(24): UVRC.asc(1).MVinput     	MV tracking       					[]
	0;...		%x(25): UVRC.asc(2).MVinput     	MV tracking       					[]
	1600;...	%x(26): V_001.P 					Pressure 							[kPa]
	270;...		%x(27): V_002.T  					Temperature                         [K]
	0];			%x(28): Vg,flare       				Totalized flaring gas std. vol.     [SMm3]

fUB = [...
	0.5; ...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0.5; ...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	3; ...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	3; ...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0.5; ...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	2500; ...	%x(6):  SG_001.P        			Pressure                            [kPa]
	320; ...	%x(7):  SG_001.T        			Temperature                         [K]
	3; ...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0.5; ...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	900; ...	%x(10): V_TO_001.P      			Pressure                            [kPa]
	370; ...	%x(11): V_TO_001.T      			Temperature                         [K]
	3; ...		%x(12): TO_001.hW      				Water phase level                   [m]
	370; ...	%x(13): TO_001.T       				Temperature                         [K]
	3; ...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0.1; ...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	300; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	3; ...		%x(18): TO_002.hW      				Water phase level                   [m]
	370; ...	%x(19): TO_002.T       				Temperature                         [K]
	300; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	330; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	900; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	330; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	1; ...		%x(24): UVRC.asc(1).MVinput     	MV tracking       					[]
	1; ...		%x(25): UVRC.asc(2).MVinput     	MV tracking       					[]
	2500; ...	%x(26): V_001.P 					Pressure 							[kPa]
	330; ...	%x(27): V_002.T  					Temperature                         [K]
	1e12]; 		%x(28): Vg,flare       				Totalized flaring gas std. vol.     [SMm3]

%Lower bounds states' trajectory inequality constraints
cLB=zeros(jmax,1);
for j=1:jmax
	for i=1:nk-1
		cLB(i+ko)=-(f(i*ne,j)-fLB(j));
	end
	ko=ko+i;
end

%Upper bounds states' trajectory inequality constraints
cUB=zeros(jmax,1);
for j=1:jmax
	for i=1:nk-1
		cUB(i+ko)=f(i*ne,j)-fUB(j);
	end
	ko=ko+i;
end

%outputting inequality constraints
c=[cLB; cUB];



end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
