classdef strW < ModelTypes
    methods
        function obj = strW(Name)
        %MODEL BRIEF
            obj.Description = 'Water stream basic model';
            obj.Name = Name;
           
        %VARIABLES
            AddVariable(obj,'Qw',flow_vol);
            AddVariable(obj,'TOG',admensional);
            AddVariable(obj,'P',pressure);
            AddVariable(obj,'T',temperature);
            AddVariable(obj,'eW',enth_mass);
            AddVariable(obj,'rhoW',dens_mass);
            AddVariable(obj,'rho_Wstd',dens_mass);
            AddVariable(obj,'d60',admensional);
            AddVariable(obj,'Tbn',temperature);
            AddVariable(obj,'ws',admensional);

            idx = GetVarIndex(obj,'Qw');        obj.Variables{idx}.Name = strcat(Name,'.Qw');       obj.Variables{idx}.Description = 'Water std. vol. flow rate';             	obj.Variables{idx}.Value = 1e4;     obj.Variables{idx}.Unit = 'Sm^3/d';
            idx = GetVarIndex(obj,'TOG');       obj.Variables{idx}.Name = strcat(Name,'.TOG');      obj.Variables{idx}.Description = 'Grease and oil content';               	obj.Variables{idx}.Value = 100;     obj.Variables{idx}.Unit = 'Sml/Sm^3';
            idx = GetVarIndex(obj,'P');         obj.Variables{idx}.Name = strcat(Name,'.P');        obj.Variables{idx}.Description = 'Pressure';                                obj.Variables{idx}.Value = 1000;    obj.Variables{idx}.Unit = 'kPa';
            idx = GetVarIndex(obj,'T');         obj.Variables{idx}.Name = strcat(Name,'.T');        obj.Variables{idx}.Description = 'Temperature';                             obj.Variables{idx}.Value = 300;     obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'eW');        obj.Variables{idx}.Name = strcat(Name,'.eW');       obj.Variables{idx}.Description = 'Water specific enthalpy';                 obj.Variables{idx}.Value = 250;     obj.Variables{idx}.Unit = 'kJ/kg';
            idx = GetVarIndex(obj,'rhoW');      obj.Variables{idx}.Name = strcat(Name,'.rhoW');     obj.Variables{idx}.Description = 'Water phase mass density';               	obj.Variables{idx}.Value = 1030;    obj.Variables{idx}.Unit = 'kg/m^3';
            idx = GetVarIndex(obj,'rho_Wstd');  obj.Variables{idx}.Name = strcat(Name,'.rho_Wstd'); obj.Variables{idx}.Description = 'Water phase mass density at std. cond.';  obj.Variables{idx}.Value = 1005;    obj.Variables{idx}.Unit = 'kg/m^3';
            idx = GetVarIndex(obj,'d60');       obj.Variables{idx}.Name = strcat(Name,'.d60');      obj.Variables{idx}.Description = 'Emulsified oil specific gravity 60 F';    obj.Variables{idx}.Value = 0.8;    	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'Tbn');       obj.Variables{idx}.Name = strcat(Name,'.Tbn');      obj.Variables{idx}.Description = 'Emulsified oil normal boiling point';     obj.Variables{idx}.Value = 650;     obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'ws');        obj.Variables{idx}.Name = strcat(Name,'.ws');       obj.Variables{idx}.Description = 'Water salt mass fraction';                obj.Variables{idx}.Value = 0.2;     obj.Variables{idx}.Unit = '';
        end
    end
end