classdef triphase_sink < ModelTypes
    methods
        function obj = triphase_sink(Name)
        %MODEL BRIEF
            obj.Description = '3-phase source stream model';
            obj.Name = Name;

        %EQUIPMENT
            AddEquipment(obj,str3ph(strcat(obj.Name,'.Inlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        end
    end
end

