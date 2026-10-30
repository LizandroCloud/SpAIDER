classdef water_source < ModelTypes
    methods
        function obj = water_source(Name)
        %MODEL BRIEF
            obj.Description = 'Water source stream model';
            obj.Name = Name;

		%Variables
            AddVariable(obj,'Cws',admensional);

            idx = GetVarIndex(obj,'Cws');    obj.Variables{idx}.Name = strcat(Name,'.Cws');   obj.Variables{idx}.Description = 'Water phase salt mass concentration';   obj.Variables{idx}.Value = 200;  obj.Variables{idx}.Unit = 'g/l';

        %EQUIPMENT
            AddEquipment(obj,water_stream(strcat(obj.Name,'.Outlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
            
        %EQUATIONS
       %Variable tokens
            Cws = strcat(obj.Name,'.Cws');
            Outlet.ws = strcat(obj.Name,'.Outlet.ws');
        
        %Water phase (brine) salt mass concentration
            AddEquation(obj,strcat(Cws,'*(1-',Outlet.ws,')',' - PP.rho_H2Ostd*',Outlet.ws));
        end
    end
end

