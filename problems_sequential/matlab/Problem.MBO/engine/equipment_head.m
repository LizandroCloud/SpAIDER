classdef equipment_head < handle
    properties
        Equipment = {};
        NumOfEqps = 0;
    end
    
    methods
        function [] = AddEquipment(obj,EqpType)
            SetEquipment(obj,EqpType);
        end
        
        function [] = SetEquipment(obj,EqpType)
            k = GetNumOfEqps(obj);
            SetNumOfEqps(obj,k+1);
            obj.Equipment{k+1,1} = EqpType;
        end
        
        function obj = SetNumOfEqps(obj,value)
           obj.NumOfEqps = value;
        end
        
        function value = GetNumOfEqps(obj)
           value = obj.NumOfEqps;
        end
        
        function idx = GetEqpIndex(obj,EqpName)
            torf = true;
            for k=1:obj.NumOfEqps
                if strcmp(obj.Equipment{k}.Name,EqpName) == true && torf
                    idx = k;
                    torf = false;
                end
            end
        end 
        
    end

end

