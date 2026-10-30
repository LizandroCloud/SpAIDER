classdef specifications_head < handle
    properties
        Specifications = {};
        NumOfSpecs = 0;
    end
    
    methods
        function [] = AddSpecification(obj,expr)
            SetSpecification(obj,expr);
        end
        
        function [] = SetSpecification(obj,expr)
            k = GetNumOfSpecs(obj);
            SetNumOfSpecs(obj,k+1);
            obj.Specifications{k+1,1} = expr;
        end
        
        function obj = SetNumOfSpecs(obj,value)
           obj.NumOfSpecs = value;
        end
        
        function value = GetNumOfSpecs(obj)
           value = obj.NumOfSpecs;
        end
    end

end

