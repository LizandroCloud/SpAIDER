classdef heat_rate < VariableTypes    
    methods
        function obj = heat_rate(~)
            obj.Value = 1e4;
            obj.Description = 'Power';
            obj.Unit = 'kW';
        end
    end    
end

