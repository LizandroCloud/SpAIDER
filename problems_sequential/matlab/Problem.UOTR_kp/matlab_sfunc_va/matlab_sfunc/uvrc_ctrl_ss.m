function [xss,uss]=uvrc_ctrl_ss(x,u,yss,par)
%uvrc_ctrl_ss - uvrc_ctrl model steady-state condition computation
%
%This function updates the values in x and u after solving the 
%steady-state condition for the uvrc_ctrl.m S-Function file.
%
%x and u are, respectively, the state and input variables passed from
%the user, which are initial guess for the steady-state computation.
%yss is the steady-state condition for the uvrc_ctrl.m outputs iden-
%tified as controlled variables after building a multi-loop control 
%scheme. par are the uvrc.m model parameters.
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
%x(8): 			TIC(1).derivFilt   		Derivative action with filter      	[]
%x(11): 		TIC(2).derivFilt   		Derivative action with filter      	[]
%x(14):  		MCC.derivFilt   		Derivative action with filter      	[]

%Free:
%x(3):  		scrb(2).P     			Pressure       						[kPa]
%x(5):          asc(1).MVinput     		MV tracking       					[]
%x(6):          asc(2).MVinput     		MV tracking       					[]
%x(7):  		TIC(1).AWBTintterm 		Integral action with AWBT			[]
%x(9): 			TIC(1).MVinput     		MV tracking       					[]
%x(10):  		TIC(2).AWBTintterm 		Integral action with AWBT			[]
%x(12): 		TIC(2).MVinput     		MV tracking       					[]
%x(13):  		MCC.AWBTintterm 		Integral action with AWBT			[]
%x(15):  		MCC.MVinput     		MV tracking       					[]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%Fixed:
%u(1):  		feed(1).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(2):  		feed(1).T     			Feed temperature       				[K]
%u(3):  		feed(2).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(4):  		feed(2).T     			Feed temperature       				[K]
%u(5):			cool(1).DP   			Pressure drop        				[kPa]
%u(6):			cool(2).DP   			Pressure drop        				[kPa]
%u(7):			P_out   				Outlet pressure        				[kPa]

%Free:

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
	false; ...	asc(1).MVinput
	false; ...	asc(2).MVinput
	false; ...	TIC(1).AWBTintterm
	true; ...	TIC(1).derivFilt
	false; ...	TIC(1).MVinput
	false; ...	TIC(2).AWBTintterm
	true; ...	TIC(2).derivFilt
	false; ...	TIC(2).MVinput
	false; ...	MCC.AWBTintterm
	true; ...	MCC.derivFilt
	false ...	MCC.MVinput
	];

fixed_x	= x(idx==true);
free_x 	= x(idx==false);

%----- Fixed and free inputs:
idu = [...
	true; ... 	feed(1).Qg
	true; ... 	feed(1).T
	true; ... 	feed(2).Qg
	true; ... 	feed(2).T
	true; ... 	cool(1).DP
	true; ... 	cool(2).DP
	true ... 	P_out
	];
fixed_u = u(idu==true);
free_u 	= u(idu==false);

%----- Fixed outputs:
idy 	= zeros(2*par.nS+2*par.nS+par.nS+2,1);

%Steady-state computation
opt = optimoptions('fsolve','Display','iter','MaxFunEvals',1e5,'MaxIter',500);
%opt = optimoptions('fsolve','Algorithm','trust-region-reflective','Display','iter','MaxFunEvals',1e5,'MaxIter',500);
z = fsolve(...
	@(z)fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par), ...
	[free_x;free_u],opt);
if imag(z)>1e-4;
	error('solution is imag')
end
free_x = real(z(1:size(find(idx==false),1)));
free_u = real(z(size(find(idx==false),1)+1:end));
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
	
	res = uvrc_ctrl([],x,u,1,[],par);
	res = res([(1:6)';9;12;15]);
	%res = res([(1:6)';[7;8];[10;11];[13;14]]);
	%res = res([(1:6)';[8;9];[11;12];[14;15]]);
	%res = res([(1:6)';7;10;13]);
    if ~isempty(yss)
        y = uvrc([],x,u,3,[],par);
        res = [res;(yss-y(idy==true))];
    end
end