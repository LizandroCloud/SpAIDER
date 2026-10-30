classdef strGL < ModelTypes
    methods
        function obj = strGL(Name)
        %MODEL BRIEF
            obj.Description = 'Gas-liquid stream basic model';
            obj.Name = Name;
           
        %VARIABLES
            AddVariable(obj,'Qg',flow_vol);
            AddVariable(obj,'Qo',flow_vol);
            AddVariable(obj,'P',pressure);
            AddVariable(obj,'T',temperature);
            AddVariable(obj,'eO',enth_mass);
            AddVariable(obj,'eG',enth_mass);
            AddVariable(obj,'rhoL',dens_mass);
            AddVariable(obj,'rhoG',dens_mass);
            AddVariable(obj,'Rs',admensional);
            AddVariable(obj,'Rv',admensional);
            AddVariable(obj,'Bo',admensional);
            AddVariable(obj,'Bg',admensional);
            AddVariable(obj,'d60',admensional);
            AddVariable(obj,'d60G',admensional);
            AddVariable(obj,'Tbn',temperature);
            AddVariable(obj,'TbnG',temperature);
            AddVariable(obj,'dg',admensional);
            AddVariable(obj,'dgs',admensional);
            AddVariable(obj,'Ycg',admensional);
            AddVariable(obj,'Ycgs',admensional);

            idx = GetVarIndex(obj,'Qg');        obj.Variables{idx}.Name = strcat(Name,'.Qg');    obj.Variables{idx}.Description = 'Gas std. vol. flow rate';               obj.Variables{idx}.Value = 1e6;     obj.Variables{idx}.Unit = 'Sm^3/d';
            idx = GetVarIndex(obj,'Qo');        obj.Variables{idx}.Name = strcat(Name,'.Qo');    obj.Variables{idx}.Description = 'Oil std. vol. flow rate';               obj.Variables{idx}.Value = 1e4;     obj.Variables{idx}.Unit = 'Sm^3/d';
            idx = GetVarIndex(obj,'P');         obj.Variables{idx}.Name = strcat(Name,'.P');     obj.Variables{idx}.Description = 'Pressure';                              obj.Variables{idx}.Value = 1000;    obj.Variables{idx}.Unit = 'kPa';
            idx = GetVarIndex(obj,'T');         obj.Variables{idx}.Name = strcat(Name,'.T');     obj.Variables{idx}.Description = 'Temperature';                           obj.Variables{idx}.Value = 300;     obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'eO');        obj.Variables{idx}.Name = strcat(Name,'.eO');    obj.Variables{idx}.Description = 'Oil specific enthalpy';                 obj.Variables{idx}.Value = 250;     obj.Variables{idx}.Unit = 'kJ/kg';
            idx = GetVarIndex(obj,'eG');        obj.Variables{idx}.Name = strcat(Name,'.eG');    obj.Variables{idx}.Description = 'Gas specific enthalpy';                 obj.Variables{idx}.Value = 400;     obj.Variables{idx}.Unit = 'kJ/kg';
            idx = GetVarIndex(obj,'rhoL');      obj.Variables{idx}.Name = strcat(Name,'.rhoL');  obj.Variables{idx}.Description = 'Liquid phase mass density';             obj.Variables{idx}.Value = 900;     obj.Variables{idx}.Unit = 'kg/m^3';
            idx = GetVarIndex(obj,'rhoG');      obj.Variables{idx}.Name = strcat(Name,'.rhoG');  obj.Variables{idx}.Description = 'Gas phase mass density';                obj.Variables{idx}.Value = 10;      obj.Variables{idx}.Unit = 'kg/m^3';
            idx = GetVarIndex(obj,'Rs');        obj.Variables{idx}.Name = strcat(Name,'.Rs');    obj.Variables{idx}.Description = 'Solution gas-oil ratio';                obj.Variables{idx}.Value = 10;      obj.Variables{idx}.Unit = 'Sm^3/Sm^3';
            idx = GetVarIndex(obj,'Rv');        obj.Variables{idx}.Name = strcat(Name,'.Rv');    obj.Variables{idx}.Description = 'Vaporized oil-gas ratio';               obj.Variables{idx}.Value = 0.01;    obj.Variables{idx}.Unit = 'Sm^3/SMm^3';
            idx = GetVarIndex(obj,'Bo');        obj.Variables{idx}.Name = strcat(Name,'.Bo');    obj.Variables{idx}.Description = 'Oil phase formation volume factor';     obj.Variables{idx}.Value = 1.05;    obj.Variables{idx}.Unit = 'm^3/Sm^3';
            idx = GetVarIndex(obj,'Bg');        obj.Variables{idx}.Name = strcat(Name,'.Bg');    obj.Variables{idx}.Description = 'Gas phase formation volume factor';     obj.Variables{idx}.Value = 0.1;     obj.Variables{idx}.Unit = 'm^3/Sm^3';
            idx = GetVarIndex(obj,'d60');       obj.Variables{idx}.Name = strcat(Name,'.d60');   obj.Variables{idx}.Description = 'Oil specific gravity 60 F';             obj.Variables{idx}.Value = 0.85;    obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'d60G');      obj.Variables{idx}.Name = strcat(Name,'.d60G');  obj.Variables{idx}.Description = 'Hypothetical d60 for gas phase';        obj.Variables{idx}.Value = 0.6;     obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'Tbn');       obj.Variables{idx}.Name = strcat(Name,'.Tbn');   obj.Variables{idx}.Description = 'Oil normal boiling point';              obj.Variables{idx}.Value = 450;     obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'TbnG');      obj.Variables{idx}.Name = strcat(Name,'.TbnG');  obj.Variables{idx}.Description = 'Hypothetical Tbn for gas phase';        obj.Variables{idx}.Value = 250;     obj.Variables{idx}.Unit = 'K';
            idx = GetVarIndex(obj,'dg');        obj.Variables{idx}.Name = strcat(Name,'.dg');    obj.Variables{idx}.Description = 'Gas specific gravity';                  obj.Variables{idx}.Value = 0.7;     obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'dgs');       obj.Variables{idx}.Name = strcat(Name,'.dgs');   obj.Variables{idx}.Description = 'Solution gas specific gravity';         obj.Variables{idx}.Value = 0.95;    obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'Ycg');       obj.Variables{idx}.Name = strcat(Name,'.Ycg');   obj.Variables{idx}.Description = 'Condensate gas-yield';                  obj.Variables{idx}.Value = 4;       obj.Variables{idx}.Unit = 'Sm^3/SMm^3';
            idx = GetVarIndex(obj,'Ycgs');      obj.Variables{idx}.Name = strcat(Name,'.Ycgs');  obj.Variables{idx}.Description = 'Soluction condensate gas-yield';        obj.Variables{idx}.Value = 5.5;     obj.Variables{idx}.Unit = 'Sm^3/SMm^3';
        end
    end
end

