classdef lowrank
    
    properties
        U
        V
    end
    
    methods
        function obj = lowrank(U, V)
            obj.U = U;
            obj.V = V;
        end
        
        function Y = mtimes(obj, X)
            Y = obj.U * (obj.V' * X);
        end
    end
end

