classdef time_sec < VariableTypes    
    methods
        function obj = time_sec(~)
            obj.Value = 60;
            obj.Description = 'Seconds';
            obj.Unit = 's';
        end
    end    
end

