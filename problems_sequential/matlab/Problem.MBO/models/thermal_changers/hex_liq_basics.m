classdef hex_liq_basics < ModelTypes  
    methods
        function obj = hex_liq_basics(Name)
        %MODEL BRIEF
            obj.Description = 'Basic equations for stream heat exchange calculations';
            obj.Name = Name;
            
        %EQUIPMENT
            AddEquipment(obj,str3ph(strcat(obj.Name,'.Inlet')));
            AddEquipment(obj,triphase_stream(strcat(obj.Name,'.Outlet')));

        %VARIABLES
            AddVariable(obj,'DP',pressure);
            AddVariable(obj,'Rseq',admensional);
            AddVariable(obj,'SumInletGas',flow_vol);
            AddVariable(obj,'Gunbal',flow_vol);

            idx = GetVarIndex(obj,'DP');    		obj.Variables{idx}.Name = strcat(Name,'.DP');   		obj.Variables{idx}.Description = 'Pressure drop';   						obj.Variables{idx}.Value = 10;  obj.Variables{idx}.Unit = 'kPa';
            idx = GetVarIndex(obj,'Rseq');    		obj.Variables{idx}.Name = strcat(Name,'.Rseq');   		obj.Variables{idx}.Description = 'Solution gas-oil ratio at eq.';      		obj.Variables{idx}.Value = 10;  obj.Variables{idx}.Unit = 'Sm^3/Sm^3';
            idx = GetVarIndex(obj,'SumInletGas');   obj.Variables{idx}.Name = strcat(Name,'.SumInletGas');  obj.Variables{idx}.Description = 'Inlet streams gas std. vol. flow sum';   	obj.Variables{idx}.Value = 2e5;	obj.Variables{idx}.Unit = 'Sm^3/d';
            idx = GetVarIndex(obj,'Gunbal');    	obj.Variables{idx}.Name = strcat(Name,'.Gunbal');   	obj.Variables{idx}.Description = 'Gas std. vol. unbalance'; 				obj.Variables{idx}.Value = 2e5; obj.Variables{idx}.Unit = 'Sm^3/d';

        %Variable tokens
            DP          = strcat(obj.Name,'.DP');
            Rseq        = strcat(obj.Name,'.Rseq');
            SumInletGas	= strcat(obj.Name,'.SumInletGas');
            Gunbal      = strcat(obj.Name,'.Gunbal');

            Inlet.Qf            = strcat(obj.Name,'.Inlet.Qf');
            Inlet.Qg            = strcat(obj.Name,'.Inlet.Qg'); 
            Inlet.BSW           = strcat(obj.Name,'.Inlet.BSW');
            Inlet.P             = strcat(obj.Name,'.Inlet.P');
            Inlet.T             = strcat(obj.Name,'.Inlet.T');
            Inlet.eO            = strcat(obj.Name,'.Inlet.eO');
            Inlet.eG            = strcat(obj.Name,'.Inlet.eG');
            Inlet.eW            = strcat(obj.Name,'.Inlet.eW');
            Inlet.rho_Wstd      = strcat(obj.Name,'.Inlet.rho_Wstd');
            Inlet.Rs            = strcat(obj.Name,'.Inlet.Rs');
            Inlet.Rv            = strcat(obj.Name,'.Inlet.Rv');
            Inlet.d60           = strcat(obj.Name,'.Inlet.d60');
            Inlet.d60G          = strcat(obj.Name,'.Inlet.d60G');
            Inlet.Tbn           = strcat(obj.Name,'.Inlet.Tbn');
            Inlet.TbnG          = strcat(obj.Name,'.Inlet.TbnG');
            Inlet.dg            = strcat(obj.Name,'.Inlet.dg');
            Inlet.dgs           = strcat(obj.Name,'.Inlet.dgs');
            Inlet.Ycg           = strcat(obj.Name,'.Inlet.Ycg');
            Inlet.Ycgs          = strcat(obj.Name,'.Inlet.Ycgs');
            Inlet.ws            = strcat(obj.Name,'.Inlet.ws');

            Outlet.Qf            = strcat(obj.Name,'.Outlet.Qf');
            Outlet.Qg            = strcat(obj.Name,'.Outlet.Qg'); 
            Outlet.BSW           = strcat(obj.Name,'.Outlet.BSW');
            Outlet.P             = strcat(obj.Name,'.Outlet.P');
            Outlet.T             = strcat(obj.Name,'.Outlet.T');
            Outlet.eO            = strcat(obj.Name,'.Outlet.eO');
            Outlet.eG            = strcat(obj.Name,'.Outlet.eG');
            Outlet.eW            = strcat(obj.Name,'.Outlet.eW');
            Outlet.rho_Wstd      = strcat(obj.Name,'.Outlet.rho_Wstd');
            Outlet.Rs            = strcat(obj.Name,'.Outlet.Rs');
            Outlet.Rv            = strcat(obj.Name,'.Outlet.Rv');
            Outlet.d60           = strcat(obj.Name,'.Outlet.d60');
            Outlet.d60G          = strcat(obj.Name,'.Outlet.d60G');
            Outlet.Tbn           = strcat(obj.Name,'.Outlet.Tbn');
            Outlet.TbnG          = strcat(obj.Name,'.Outlet.TbnG');
            Outlet.dg            = strcat(obj.Name,'.Outlet.dg');
            Outlet.dgs           = strcat(obj.Name,'.Outlet.dgs');
            Outlet.Ycg           = strcat(obj.Name,'.Outlet.Ycg');
            Outlet.Ycgs          = strcat(obj.Name,'.Outlet.Ycgs');
            Outlet.ws            = strcat(obj.Name,'.Outlet.ws');     
            
        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        
        %EQUATIONS
        %Oil std. vol. balance
            AddEquation(obj,strcat(...
                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,') + ',Inlet.Qg,'/(24*3600)*',Inlet.Rv,'*1e-3',' - (', ...
					Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,') + ',Outlet.Qg,'/(24*3600)*',Outlet.Rv,'*1e-3', ...
                ')'));

        %Sum of inlet streams gas std. vol.
            AddEquation(obj,strcat( ...
				SumInletGas,'/(24*3600)',' - (',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Rs,' + ',Inlet.Qg,'/(24*3600)',')'));

        %Gas std. vol. unbalance
            AddEquation(obj,strcat(...
                Gunbal,'/(24*3600)',' - (',Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Rseq,' - ',SumInletGas,'/(24*3600)',')'));

		%Std. vol. balance for liberated gas
            AddEquation(obj,strcat(...
                Outlet.Qg,'/(24*3600)',' - 0.5*(1+tanh(-',Gunbal,'/1e-6))*(-',Gunbal,'/(24*3600)',')'));

        %Std. vol. balance for solution gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Outlet.Rs,' - (', ...
					'0.5*(1+tanh(',Gunbal,'/1e-6))*',SumInletGas,'/(24*3600)', ...
					' + 0.5*(1+tanh(-',Gunbal,'/1e-6))*',Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Rseq, ...
				')'));
            
        %Water std. vol. balance
            AddEquation(obj,strcat(Inlet.Qf,'/(24*3600)*',Inlet.BSW,' - ',Outlet.Qf,'/(24*3600)*',Outlet.BSW));
        
        %Water salt mass balance
            AddEquation(obj,strcat(Inlet.ws,' - ',Outlet.ws));

        %Pressure drop
            AddEquation(obj,strcat(DP,' - (',Inlet.P,' - ',Outlet.P,')'));
            
        %Solution gas-oil ratio at equilibrium
            AddEquation(obj,strcat(...
                Rseq,' - SolGasOilRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,')'));


        %Vaporized gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Rv,' - VapOilGasRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,',',Outlet.Ycg,')'));

        %Equates of outlet stream compositional variables
            AddEquation(obj,strcat(Outlet.d60,' - ',Inlet.d60));
            AddEquation(obj,strcat(Outlet.Tbn,' - ',Inlet.Tbn));
			AddEquation(obj,strcat(Outlet.dg,' - ',Inlet.dg));
			AddEquation(obj,strcat(Outlet.Ycg,' - ',Inlet.Ycg));

        end
    end
    
end

