classdef gasliq_sink < ModelTypes
    methods
        function obj = gasliq_sink(Name)
        %MODEL BRIEF
            obj.Description = '3-phase source stream model';
            obj.Name = Name;

        %EQUIPMENT
            AddEquipment(obj,strGL(strcat(obj.Name,'.Inlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        end
    end
end

