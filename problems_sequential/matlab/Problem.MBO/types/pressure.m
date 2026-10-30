classdef pressure < VariableTypes    
    methods
        function obj = pressure(~)
            obj.Value = 100;
            obj.Description = 'Pressure';
            obj.Unit = 'kPa';
        end
    end    
end

