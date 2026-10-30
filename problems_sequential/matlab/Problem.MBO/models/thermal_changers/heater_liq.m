classdef heater_liq < hex_liq_basics  
    methods
        function obj = heater_liq(Name)
        %MODEL BRIEF
            obj = obj@hex_liq_basics(Name);
			obj.Description = 'Basic equations for stream heating calculations';
            obj.Name = Name;

        %VARIABLES
            AddVariable(obj,'Qin',heat_rate);

            idx = GetVarIndex(obj,'Qin');	obj.Variables{idx}.Name = strcat(Name,'.Qin');	obj.Variables{idx}.Description = 'Heat applied to the system';	obj.Variables{idx}.Value = 16e3;	obj.Variables{idx}.Unit = 'kW';

        %Variable tokens
            Qin          	= strcat(obj.Name,'.Qin');

            Inlet.Qf      	= strcat(obj.Name,'.Inlet.Qf');
            Inlet.Qg      	= strcat(obj.Name,'.Inlet.Qg'); 
            Inlet.BSW     	= strcat(obj.Name,'.Inlet.BSW');
            Inlet.P       	= strcat(obj.Name,'.Inlet.P');
            Inlet.T       	= strcat(obj.Name,'.Inlet.T');
            Inlet.eO      	= strcat(obj.Name,'.Inlet.eO');
            Inlet.eG      	= strcat(obj.Name,'.Inlet.eG');
            Inlet.eW      	= strcat(obj.Name,'.Inlet.eW');
            Inlet.rho_Wstd	= strcat(obj.Name,'.Inlet.rho_Wstd');
            Inlet.Rs      	= strcat(obj.Name,'.Inlet.Rs');
            Inlet.Rv      	= strcat(obj.Name,'.Inlet.Rv');
            Inlet.d60     	= strcat(obj.Name,'.Inlet.d60');
            Inlet.d60G    	= strcat(obj.Name,'.Inlet.d60G');
            Inlet.Tbn     	= strcat(obj.Name,'.Inlet.Tbn');
            Inlet.TbnG    	= strcat(obj.Name,'.Inlet.TbnG');
            Inlet.dg      	= strcat(obj.Name,'.Inlet.dg');
            Inlet.dgs     	= strcat(obj.Name,'.Inlet.dgs');
            Inlet.Ycg     	= strcat(obj.Name,'.Inlet.Ycg');
            Inlet.Ycgs    	= strcat(obj.Name,'.Inlet.Ycgs');
            Inlet.ws      	= strcat(obj.Name,'.Inlet.ws');

            Outlet.Qf      	= strcat(obj.Name,'.Outlet.Qf');
            Outlet.Qg      	= strcat(obj.Name,'.Outlet.Qg'); 
            Outlet.BSW     	= strcat(obj.Name,'.Outlet.BSW');
            Outlet.P       	= strcat(obj.Name,'.Outlet.P');
            Outlet.T       	= strcat(obj.Name,'.Outlet.T');
            Outlet.eO      	= strcat(obj.Name,'.Outlet.eO');
            Outlet.eG      	= strcat(obj.Name,'.Outlet.eG');
            Outlet.eW      	= strcat(obj.Name,'.Outlet.eW');
            Outlet.rho_Wstd	= strcat(obj.Name,'.Outlet.rho_Wstd');
            Outlet.Rs      	= strcat(obj.Name,'.Outlet.Rs');
            Outlet.Rv      	= strcat(obj.Name,'.Outlet.Rv');
            Outlet.d60     	= strcat(obj.Name,'.Outlet.d60');
            Outlet.d60G    	= strcat(obj.Name,'.Outlet.d60G');
            Outlet.Tbn     	= strcat(obj.Name,'.Outlet.Tbn');
            Outlet.TbnG    	= strcat(obj.Name,'.Outlet.TbnG');
            Outlet.dg      	= strcat(obj.Name,'.Outlet.dg');
            Outlet.dgs     	= strcat(obj.Name,'.Outlet.dgs');
            Outlet.Ycg     	= strcat(obj.Name,'.Outlet.Ycg');
            Outlet.Ycgs    	= strcat(obj.Name,'.Outlet.Ycgs');
            Outlet.ws      	= strcat(obj.Name,'.Outlet.ws');     
            
        %REVERSE MODEL INHERITANCE
%            SystemVariables(obj);
%            SystemEquations(obj);
%            SystemInitConditions(obj);
        
        %EQUATIONS
        %Energy balance
        %Energy balance
            AddEquation(obj,strcat(...
                Qin,' + ',Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd*',Inlet.eO, ...
						' + ',Inlet.BSW,'*',Inlet.rho_Wstd,'*',Inlet.eW,')', ...
                    ' + ',Inlet.Qg,'/(24*3600)*',Inlet.dg,'*PP.rho_AIRstd*',Inlet.eG,' - (', ...
                Outlet.Qf,'/(24*3600)*((1-',Outlet.BSW,')*',Outlet.d60,'*PP.rho_H2Ostd*',Outlet.eO, ...
						' + ',Outlet.BSW,'*',Outlet.rho_Wstd,'*',Outlet.eW,')', ...
                    ' + ',Outlet.Qg,'/(24*3600)*',Outlet.dg,'*PP.rho_AIRstd*',Outlet.eG, ...
                ')'));

        end
    end
    
end

