% function [ceq c dceq dc]=cr2(f,t_ult,t_out,u,y,stages,ns_a)
function [ceq c gradceq gradc]=cr2(f,t_ult,t_out,df,u,y,stages,ns_a)
global modelo_u0 idu
% ts = stages;
ts =mat2cell(stages, 1,length(stages));
% for 1:length(u)
    
%jmax= length(f(1,:));       	%states
%nk	= length(stages)-1;     	%number of elements
%ne	= (length(f(:,1))-1)/nk;	%length of the element
Ndec= size(u(:),1);				%Number of decision variables
Nt	= size(t_out,1);			%Number of time steps
n_c = length(ns_a);
ws = cell(n_c,1);
% ws{j} = mat2cell(u, n_c, length(u));
k = 1;
for j=1:n_c
    ws{j} = u(k:k+ns_a(j)-1);
    k = k+ns_a(j);
end
%jmax= length(f(1,:));       	%states
%nk	= length(stages)-1;     	%number of elements
%ne	= (length(f(:,1))-1)/nk;	%length of the element
% Ndec= size(cat(1,ws{:}),1);  %Number of decision variables
% n_c = size(ts,2);               %Number of control variables
% Nt	= size(t_out,1);            %Number of time steps

ceq=[]; c=[]; gradceq=[]; gradc=[];

%----- Equality constraints ---------------------------------------
%Equality constraint for u initial condition
j0=0;
for j=1:n_c
    if j==1
        ceqi(j)=ws{j}(1)-modelo_u0(idu(j))*1e-3;
    else
        ceqi(j)=ws{j}(1)-modelo_u0(idu(j));
    end
    gradceqi(:,j)=zeros(Ndec,1);
    gradceqi(j0+1,j)=1;
    j0=j0+length(ws{j});
end

%Equality constraint for u(1) for horizon control
idts=min(find(ts{1}>=1800));
ceqh=ws{1}(idts:end)-0.2*modelo_u0(idu(1))*1e-3;
gradceqi(idts:length(ws{1}),length(ceqi)+1:length(ceqi)+length(ceqh))=eye(length(ceqh));

% k = 1;
% for j=1:size(ns_a,1)
%     if j==1
%         ceqi(j)=u(k)-modelo_u0(idu(j))*1e-3;
%     else
%         ceqi(j)=u(k)-modelo_u0(idu(j));
%     end
%     k = k+ns_a(j);
% end

%outputting equality constraints
% ceq=ceqi(:);
% gradceq=gradceqi;
ceq=[ceqi(:); ceqh(:)];
gradceq=gradceqi;
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
%cLB=-(f(:,idx)-ones(Nt,1)*fLB'); cLB=cLB(:);
cLB=-(f(:,idx)-ones(Nt,1)*fLB')./(ones(Nt,1)*fLB'); cLB=cLB(:);
%cLB=-(f(:,idx)-ones(Nt,1)*fLB(idx)'); cLB=cLB(:);

gradcLB=[];
for j=1:length(idx)
    g=[];
	for m=1:Ndec
        %Id=-df(:,idx(j),m);
        Id=-df(:,idx(j),m)./fLB(j);
		g=[g;Id(:)'];
	end
	gradcLB=[gradcLB,g];
end
%}
%Upper bounds states' trajectory inequality constraints
%cUB=(f(:,idx)-ones(Nt,1)*fUB'); cUB=cUB(:);
cUB=(f(:,idx)-ones(Nt,1)*fUB')./(ones(Nt,1)*fUB'); cUB=cUB(:);
%cUB=(f(:,idx)-ones(Nt,1)*fUB(idx)'); cUB=cUB(:);
%gradcUB=-gradcLB;

gradcUB=[];
for j=1:length(idx)
    g=[];
	for m=1:Ndec
        %Id=-df(:,idx(j),m);
        Id=df(:,idx(j),m)./fUB(j);
		g=[g;Id(:)'];
	end
	gradcUB=[gradcUB,g];
end

%{
%Lower bounds states' trajectory inequality constraints
cLB=[]; ko=numel(cLB);
for j=1:length(idx)
    for i=1:nk-1
		cLB(i+ko)=-(f(i*ne,idx(j))-fLB(j));
        Id=-df(i*ne,idx(j),:);
        gradcLB(:,i+ko)=Id(:);
    end
	ko=ko+i;
end


%Upper bounds states' trajectory inequality constraints
cUB=[]; ko=numel(cUB);
for j=1:length(idx)
	for i=1:nk-1
		cUB(i+ko)=f(i*ne,idx(j))-fUB(j);
        Id=df(i*ne,idx(j),:);
        gradcUB(:,i+ko)=Id(:);        
	end
	ko=ko+i;
end
%}

%outputting inequality constraints
c=[cLB; cUB];
gradc=[gradcLB, gradcUB];


%Lower bounds inputs' change rate trajectory inequality constraints
dwsLB=zeros(length(ws{1})-1,1);
graddwsLB=zeros(length(ws{1})-1,1);
for i=1:length(dwsLB)
    dwsLB(i)=( (ws{1}(i)-ws{1}(i+1)) - 0.1*modelo_u0(1)*1e-3*(ts{1}(i+1)-ts{1}(i))/120 );
    graddwsLB(i+1,i)=-1;
    graddwsLB(i,i)=1;
end
dwsLB=dwsLB(:);
graddwsLB(length(dwsLB)+1:Ndec,:)=0;

%Upper bounds inputs' change rate trajectory inequality constraints
dwsUB=zeros(length(ws{1})-1,1);
graddwsUB=zeros(length(ws{1})-1,1);
for i=1:length(dwsUB)
    dwsUB(i)=( (ws{1}(i+1)-ws{1}(i)) - 0.1*modelo_u0(1)*1e-3*(ts{1}(i+1)-ts{1}(i))/120 );
    graddwsUB(i+1,i)=1;
    graddwsUB(i,i)=-1;
end
dwsUB=dwsUB(:);
graddwsUB(length(dwsUB)+1:Ndec,:)=0;

%outputting inequality constraints and gradients
c=[c; dwsLB; dwsUB];
gradc=[gradc, graddwsLB, graddwsUB];
%}
end