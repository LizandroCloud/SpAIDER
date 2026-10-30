% function [ceq,c] = cr1(f,t_ult,t_out,ts,df,u,y,stages)
function [ceq,c] = cr1(f,t_ult,t_out,u,y,stages, ns_a)
global modelo_u0
global ufinal

idu=[1;16;17;18];

ts = stages;
%jmax= length(f(1,:));       	%states
%nk	= length(stages)-1;     	%number of elements
%ne	= (length(f(:,1))-1)/nk;	%length of the element
% ws=u;
% Ndec= size(cat(1,ws{:}),1);  %Number of decision variables
n_c = size(ts,2);               %Number of control variables
Nt	= size(t_out,1);            %Number of time steps

% ns_a=[];
% for j=1:size(stages,2)
	% ns_a(j,1)=length(stages(:,j));
% % end
% ws=mat2cell(u,ns_a,1);

ceq=[]; c=[];

%----- Equality constraints ---------------------------------------
%Equality constraint for u(1) initial 
%Equality constraints for initial conditions
k = 1;
for j=1:size(ns_a,1)
    if j==1
        ceqi(j)=u(k)-modelo_u0(idu(j))*1e-3;
    else
        ceqi(j)=u(k)-modelo_u0(idu(j));
    end
    k = k+ns_a(j);
end

%outputting equality constraints
ceq=ceqi(:);
%--------------------------------------------------------------------

%----- Inequality constraints ---------------------------------------
%States upper and lower bounds
fLB = [...
	...0.1;...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	...0;...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	...1;...	%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1;...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	...0;...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	1600;...	%x(6):  SG_001.P        			Pressure                            [kPa]
	...270;...	%x(7):  SG_001.T        			Temperature                         [K]
	0.5;...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	...0;...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	700;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	...340;...	%x(11): V_TO_001.T      			Temperature                         [K]
	...0.1;...	%x(12): TO_001.hW      				Water phase level                   [m]
	...340;...	%x(13): TO_001.T       				Temperature                         [K]
	0.75;...0.5;...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	...0;...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	150;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	...340;...	%x(17): V_TO_002.T      			Temperature                         [K]
	...0.1;...	%x(18): TO_002.hW      				Water phase level                   [m]
	...340;...	%x(19): TO_002.T       				Temperature                         [K]
	150;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	...290;...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	700;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	...290;...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	1600];...	%x(27): V_001.P 					Pressure 							[kPa]
	...270];    %x(28): V_001.T  					Temperature                         [K]

fUB = [...
	...0.5; ...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	...0.5; ...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	...3; ...	%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	3; ...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	...0.5; ...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	2500; ...	%x(6):  SG_001.P        			Pressure                            [kPa]
	...320; ...	%x(7):  SG_001.T        			Temperature                         [K]
	3; ...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	...0.5; ...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	900; ...	%x(10): V_TO_001.P      			Pressure                            [kPa]
	...370; ...	%x(11): V_TO_001.T      			Temperature                         [K]
	...3; ...	%x(12): TO_001.hW      				Water phase level                   [m]
	...370; ...	%x(13): TO_001.T       				Temperature                         [K]
	1.25;...3; ...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	...0.1; ...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	320; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	...370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	...3; ...	%x(18): TO_002.hW      				Water phase level                   [m]
	...370; ...	%x(19): TO_002.T       				Temperature                         [K]
	320; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	...330; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	900; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	...330; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	2500]; ...	%x(27): V_001.P 					Pressure 							[kPa]
	...330];    %x(28): V_001.T  					Temperature                         [K]

%idx=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;27;28];
idx=[4;6;8;10;14;16;20;22;27];
%idx=4;

%Lower bounds states' trajectory inequality constraints
cLB=-(f(:,idx)-ones(Nt,1)*fLB'); cLB=cLB(:);
%cLB=-(f(:,idx)-ones(Nt,1)*fLB(idx)'); cLB=cLB(:);

%Upper bounds states' trajectory inequality constraints
cUB=(f(:,idx)-ones(Nt,1)*fUB'); cUB=cUB(:);
%cUB=(f(:,idx)-ones(Nt,1)*fUB(idx)'); cUB=cUB(:);

%{
%Lower bounds states' trajectory inequality constraints
cLB=[]; ko=numel(cLB);
for j=1:length(idx)
    for i=1:nk-1
		cLB(i+ko)=-(f(i*ne,idx(j))-fLB(j));
    end
	ko=ko+i;
end


%Upper bounds states' trajectory inequality constraints
cUB=[]; ko=numel(cUB);
for j=1:length(idx)
	for i=1:nk-1
		cUB(i+ko)=f(i*ne,idx(j))-fUB(j);
	end
	ko=ko+i;
end
%}

%outputting inequality constraints
c=[cLB; cUB];

%Lower bounds inputs' change rate trajectory inequality constraints
% njh = sum(ns_a);
dwsLB=zeros(ns_a(1)-1,1);
k=1;
    for i=1:ns_a(1)-1
        dwsLB(i)=( ( u(k)-u(k+1) ) - 0.1*modelo_u0(1)*1e-3*(ts(k+1)-ts(k))/120 );
        k=k+1;
    end

dwsLB=dwsLB(:);

%Upper bounds inputs' change rate trajectory inequality constraints
dwsUB=zeros(ns_a(1)-1,1);
k=1;
    for i=1:ns_a(1)-1
        dwsUB(i)=( ( u(k+1)-u(k) ) - 0.1*modelo_u0(1)*1e-3*(ts(k+1)-ts(k))/120 );
        k = k+1;
    end


% dwsUB=zeros(length(ws{1})-1,1);
% for i=1:length(dwsUB)
%     dwsUB(i)=( (ws{1}(i+1)-ws{1}(i)) - 0.1*modelo_u0(1)*1e-3*(ts{1}(i+1)-ts{1}(i))/120 );
% end
dwsUB=dwsUB(:);

%outputting inequality constraints and gradients
c=[c; dwsLB; dwsUB];

end