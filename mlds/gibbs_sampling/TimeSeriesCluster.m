classdef TimeSeriesCluster < handle
    %TIMESERIESCLUSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Ids
        X
        Ts
        V
        Delta_t
        Dim % dimension of the data
        
        % suffcient statistics:
        %   XX1 = \sum(x*x')
        %   VX1 = \sum(v*x')
        %   V1V = \sum(v'*v)
        %   TiMinus1 = \sum(Ti-1)
        SuffStats = []
        
        Params
    end
    
    methods
        function obj = TimeSeriesCluster(ids, X, ts, dim, params)
            obj.Ids = ids;
            obj.X = X;
            obj.Ts = ts;
            obj.Dim = dim;
            obj.Params = params;
            
            [V, delta_t] = calculate_velocity(X, ts);
            obj.V = V;
            obj.Delta_t = delta_t;
            
            suffStats = obj.calcSufficientStats(X, V);
            obj.SuffStats = suffStats;        
        end
        
        function suffStats = calcSufficientStats(obj, X, V)
            D = obj.Dim;
            XX1_new = zeros(D, D);
            VX1_new = zeros(D, D);
            V1V_new = 0;
            TiMinus1_new = 0;
            for i = 1:length(X)
                x1 = X{i}(:,1:end-1);
                XX1_new = XX1_new + x1 * x1';
                VX1_new = VX1_new + V{i} * x1';
                V1V_new = V1V_new + sum(sum(V{i}.^2));
                TiMinus1_new = TiMinus1_new + size(V{i},2);
            end
            
            suffStats = [];
            suffStats.XX1 = XX1_new;
            suffStats.VX1 = VX1_new;
            suffStats.V1V = V1V_new;
            suffStats.TiMinus1 = TiMinus1_new;
        end
        
        function addData(obj, id, x, t)
            % x and t are a cell array of multiple data
            obj.Ids = [obj.Ids, id];
            obj.X = [obj.X, x];
            obj.Ts = [obj.Ts, t];
            
            [v, delta_t] = calculate_velocity(x, t);
            obj.V = [obj.V, v];
            obj.Delta_t = [obj.Delta_t, delta_t];
            
            % add sufficient statistics
            suffStats = obj.calcSufficientStats(x, v);
            obj.addSuffStats(suffStats);
        end
        
        function addSuffStats(obj, suffStats)
            obj.SuffStats.XX1 = obj.SuffStats.XX1 + suffStats.XX1;
            obj.SuffStats.VX1 = obj.SuffStats.VX1 + suffStats.VX1;
            obj.SuffStats.V1V = obj.SuffStats.V1V + suffStats.V1V;
            obj.SuffStats.TiMinus1 = obj.SuffStats.TiMinus1 + suffStats.TiMinus1;
        end
        
        function removeData(obj, id)
            % id should be a single data to remove
            k = find(obj.Ids == id);
            suffStats = obj.calcSufficientStats(obj.X(k), obj.V(k));
            
            selected = (obj.Ids ~= id);
            
            obj.Ids = obj.Ids(selected);
            obj.X = obj.X(selected);
            obj.Ts = obj.Ts(selected);
            obj.V = obj.V(selected);
            obj.Delta_t = obj.Delta_t(selected);
            
            % remove sufficient stats
            obj.removeSuffStats(suffStats);
            
%             if length(obj.X) <= 1
%                 error('There is only one point left in the cluster');
%             end
        end
        
        function removeSuffStats(obj, suffStats)
            obj.SuffStats.XX1 = obj.SuffStats.XX1 - suffStats.XX1;
            obj.SuffStats.VX1 = obj.SuffStats.VX1 - suffStats.VX1;
            obj.SuffStats.V1V = obj.SuffStats.V1V - suffStats.V1V;
            obj.SuffStats.TiMinus1 = obj.SuffStats.TiMinus1 - suffStats.TiMinus1;
        end
        
        function val = calcLogPredictive(obj, x, t)
            error('Implement calcLogPredictive in the subclass');
        end
        
        function [A, sigma2] = sampleAsAndSigmas(obj)
            error('Implement sampleAsAndSigmas in the subclass');
        end
        
    end
end

