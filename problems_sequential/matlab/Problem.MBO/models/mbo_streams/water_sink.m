classdef water_sink < ModelTypes
    methods
        function obj = water_sink(Name)
        %MODEL BRIEF
            obj.Description = 'Water sink stream model';
            obj.Name = Name;

		%Variables
            AddVariable(obj,'Cws',admensional);

            idx = GetVarIndex(obj,'Cws');    obj.Variables{idx}.Name = strcat(Name,'.Cws');   obj.Variables{idx}.Description = 'Water phase salt mass concentration';   obj.Variables{idx}.Value = 200;  obj.Variables{idx}.Unit = 'g/l';

        %EQUIPMENT
            AddEquipment(obj,strW(strcat(obj.Name,'.Inlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
            
        %EQUATIONS
       %Variable tokens
            Cws = strcat(obj.Name,'.Cws');
            Inlet.ws = strcat(obj.Name,'.Inlet.ws');
        
        %Water phase (brine) salt mass concentration
            AddEquation(obj,strcat(Cws,'*(1-',Inlet.ws,')',' - PP.rho_H2Ostd*',Inlet.ws));
        end
    end
end