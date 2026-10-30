classdef equations_head < handle
    properties
        Equations = {};
        NumOfEqs = 0;
    end
    
    methods
        function [] = AddEquation(obj,expr)
            SetEquation(obj,expr);
        end
        
        function [] = SetEquation(obj,expr)
            k = GetNumOfEqs(obj);
            SetNumOfEqs(obj,k+1);
            obj.Equations{k+1,1} = expr;
        end
        
        function obj = SetNumOfEqs(obj,value)
           obj.NumOfEqs = value;
        end
        
        function value = GetNumOfEqs(obj)
           value = obj.NumOfEqs;
        end
    end

end

