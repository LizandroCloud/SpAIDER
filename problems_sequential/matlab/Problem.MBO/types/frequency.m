classdef frequency < VariableTypes    
    methods
        function obj = frequency(~)
            obj.Value = 1/60;
            obj.Description = 'Frequency';
            obj.Unit = 'Hz';
        end
    end    
end

