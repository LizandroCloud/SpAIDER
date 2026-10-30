classdef temperature < VariableTypes    
    methods
        function obj = temperature(~)
            obj.Value = 300;
            obj.Description = 'Temperature';
            obj.Unit = 'K';
        end
    end    
end

