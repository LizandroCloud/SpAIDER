classdef mass < VariableTypes    
    methods
        function obj = mass(~)
            obj.Value = 1e5;
            obj.Description = 'Mass';
            obj.Unit = 'kg';
        end
    end    
end

