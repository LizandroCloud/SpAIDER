function [xss,uss]=uvrc_ss(x,u,yss,par)
%uvrc_ss - Vapor recovery unit unit steady-state condition computation
%
%This function updates the values in x and u after solving the 
%steady-state condition for the vapor recovery unit upstream topside 
%facility model, described by the uvrc.m S-Function file.
%
%x and u are, respectively, the state and input variables passed from
%the user, which are initial guess for the steady-state computation.
%yss is the steady-state condition for the uvrc.m outputs identified
%as controlled variables after building a multi-loop control scheme.
%par are the uvrc.m model parameters.
%
%xss and uss are, respectively, the updated state and input variables
%at steady-state condition. The variables in xss updated are the ones
%which remain free after building the multi-loop control scheme. The 
%variables in uss updated are manipulated variables of the multi-loop
%control scheme.
%
%Fixed states: 	Controlled variables. Values not updated.
%Fixed inputs:	Exogenous input variables. Values not updated.
%Fixed outpus:  Controlled variables. Steady-state condition as yss.
%Free states:  	Free variables. Values updated.
%Free inputs: 	Manipulated variables. Values updated.
%
%----- States -------------------------------------------------------
%Fixed:
%x(1):  		scrb(1).P     			Pressure       						[kPa]
%x(2):  		scrb(1).T     			Temperature       					[K]
%x(4):  		scrb(2).T     			Temperature       					[K]
%x(7):          asc(1).MVinput     		MV tracking       					[]
%x(10):         asc(2).MVinput     		MV tracking       					[]

%Free:
%x(3):  		scrb(2).P     			Pressure       						[kPa]
%x(5):			asc(1).AWBTintterm     	Integral action with AWBT  			[]
%x(6):          asc(1).derivFilt     	Derivative action with filter       []
%x(8):			asc(2).AWBTintterm     	Integral action with AWBT  			[]
%x(9):          asc(2).derivFilt     	Derivative action with filter       []
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%Fixed:
%u(2):  		feed(1).T     			Feed temperature       				[K]
%u(4):  		feed(2).T     			Feed temperature       				[K]
%u(6):			cool(1).DP   			Pressure drop        				[kPa]
%u(8):			cool(2).DP   			Pressure drop        				[kPa]
%u(9):			P_out   				Outlet pressure        				[kPa]

%Free:
%u(1):  		feed(1).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(3):  		feed(2).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(5):			cool(1).Q_out   		Heat removed from the system        [kW]
%u(7):			cool(2).Q_out   		Heat removed from the system        [kW]
%u(10):			MV_in     				Input signal from controller        []
%--------------------------------------------------------------------
%------ Outputs for steady-state condition --------------------------
%Fixed:
%There are no fixed outputs
%--------------------------------------------------------------------

%2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%----- Fixed and free states:
idx = [...
	true; ...	scrb(1).P
	true; ...	scrb(1).T
	false; ...	scrb(2).P
	true; ...	scrb(2).T
	false; ...	asc(1).AWBTintterm
	false; ...	asc(1).derivFilt  
	true; ...	asc(1).MVinput    
	false; ...	asc(2).AWBTintterm
	false; ...	asc(2).derivFilt  
	true; ...	asc(2).MVinput    
	];
fixed_x	= x(idx==true);
free_x 	= x(idx==false);

%----- Fixed and free inputs:
idu = [...
	false; ... 	feed(1).Qg
	true; ... 	feed(1).T
	false; ... 	feed(2).Qg
	true; ... 	feed(2).T
	false; ... 	cool(1).Q_out
	true; ... 	cool(1).DP
	false; ... 	cool(2).Q_out
	true; ... 	cool(2).DP
	true; ... 	P_out
	false ... 	MV_in
	];
fixed_u = u(idu==true);
free_u 	= u(idu==false);

%----- Fixed outputs:
idy 	= zeros(2*par.nS+2*par.nS+par.nS+2,1);

%Steady-state computation
opt = optimoptions('fsolve','Display','iter','MaxFunEvals',1e5,'MaxIter',500);
z = fsolve(...
	@(z)fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par), ...
	[free_x;free_u],opt);
free_x = z(1:size(find(idx==false),1));
free_u = z(size(find(idx==false),1)+1:end);
xss = x; uss = u;
xss(idx==false) = free_x;
uss(idu==false) = free_u;
end

function res = fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par)
    free_x          = z(1:size(find(idx==false),1));
    free_u          = z(size(find(idx==false),1)+1:end);

    x(idx==true) 	= fixed_x;
    x(idx==false) 	= free_x;
	x = x(:);
	
    u(idu==true) 	= fixed_u;
	u(idu==false) 	= free_u;
	u = u(:);
	
	res = uvrc([],x,u,1,[],par);
    if ~isempty(yss)
        y = uvrc([],x,u,3,[],par);
        res = [res;(yss-y(idy==true))];
    end
end