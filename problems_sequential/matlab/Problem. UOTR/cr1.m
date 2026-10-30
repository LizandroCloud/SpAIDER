function [ceq,c]=cr1(f,t_ult,t_out,u,y,stages)
global modelo_u0 modelo_par

jmax= length(f(1,:));       	%states
nk	= length(stages)-1;     	%number of elements
ne	= (length(f(:,1))-1)/nk;	%length of the element

%----- Equality constraints ---------------------------------------
%Equality constraint for u(1) initial condition
ceq(1)=u(1)-modelo_u0(1)*1e-3;

%Equality constraint for u(1) terminal condition
ufinal=1e6/(modelo_u0(2)*(1-modelo_u0(3)));
ceq(2)=u(end)-ufinal*1e-3;

%outputting equality constraints
ceq=ceq(:);

%{
%For UOTR open-loop
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
ceq=[ceqdfterm; ceqyterm]; ko=numel(ceq);
for j=1:jmax
	for i=1:nk-1
		ceqNP(i+ko)=f(i*ne+1,j)-f(i*ne,j);
	end
	ko=ko+i;
end
ceqNP=ceqNP(:);

%outputting equality constraints
ceq=[ceq; ceqNP];
%}
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
	0.5;...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0;...		%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	700;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	340;...		%x(11): V_TO_001.T      			Temperature                         [K]
	0.1;...		%x(12): TO_001.hW      				Water phase level                   [m]
	340;...		%x(13): TO_001.T       				Temperature                         [K]
	0.5;...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0;...		%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	150;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	340;...		%x(17): V_TO_002.T      			Temperature                         [K]
	0.1;...		%x(18): TO_002.hW      				Water phase level                   [m]
	340;...		%x(19): TO_002.T       				Temperature                         [K]
	150;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	290;...		%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	700;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	290;...		%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	1600;...	%x(27): V_001.P 					Pressure 							[kPa]
	270];       %x(28): V_001.T  					Temperature                         [K]

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
	320; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	3; ...		%x(18): TO_002.hW      				Water phase level                   [m]
	370; ...	%x(19): TO_002.T       				Temperature                         [K]
	320; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	330; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	900; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	330; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	2500; ...	%x(27): V_001.P 					Pressure 							[kPa]
	330];       %x(28): V_001.T  					Temperature                         [K]

%{
%For UOTR open-loop
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
%}
idx=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;27;28];

%Lower bounds states' trajectory inequality constraints
cLB=[]; ko=0;
for j=1:length(idx)
	for i=1:nk-1
		cLB(i+ko)=-(f(i*ne,idx(j))-fLB(j));
	end
	ko=ko+i;
end

%Upper bounds states' trajectory inequality constraints
cUB=[]; ko=0;
for j=1:length(idx)
	for i=1:nk-1
		cUB(i+ko)=f(i*ne,idx(j))-fUB(j);
	end
	ko=ko+i;
end

%outputting inequality constraints
c=[cLB; cUB];

%{
%Lower bounds inputs' change rate trajectory inequality constraints
duLB=zeros(Nut-1,1); gradduLB=zeros(Nut-1);
for i=1:Nut-2
    duLB(i)=( (u(i)-u(i+1)) - 0.1*modelo_u0(1)*1e-3*(ut(i+1)-ut(i))/120 );
    gradduLB(i+1,i)=-1;
    gradduLB(i,i)=1;%gradduLB(i,i-1)=1;
end
duLB(Nut-1)=( (u(end)-ufinal) - 0.1*modelo_u0(1)*1e-3*(ut(end)-ut(end-1))/120 );
gradduLB(Nut-1,Nut-1)=1; %gradduLB(Nut,Nut-1)=1;
duLB=duLB(:);

%Upper bounds inputs' change rate trajectory inequality constraints
duUB=zeros(Nut-1,1); gradduUB=zeros(Nut-1);
for i=1:Nut-2
    duUB(i)=( (u(i+1)-u(i)) - 0.1*modelo_u0(1)*1e-3*(ut(i+1)-ut(i))/120 );
    gradduUB(i+1,i)=1;
    gradduUB(i,i)=-1;%gradduLB(i,i-1)=1;
end
duUB(Nut-1)=( (ufinal-u(end)) - 0.1*modelo_u0(1)*1e-3*(ut(end)-ut(end-1))/120 );
gradduUB(Nut-1,Nut-1)=-1; %gradduLB(Nut,Nut-1)=1;
duUB=duUB(:);

%outputting inequality constraints and gradients
c=[cLB; cUB; duLB; duUB];
gradc=[gradcLB, gradcUB]*1e3; %Retornando ordem de grandeza
gradc=[gradc, gradduLB, gradduUB];
%}

end