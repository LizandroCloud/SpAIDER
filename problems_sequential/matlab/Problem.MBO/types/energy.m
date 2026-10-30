classdef energy < VariableTypes    
    methods
        function obj = energy(~)
            obj.Value = 5e6;
            obj.Description = 'Energy';
            obj.Unit = 'kJ';
        end
    end    
end

