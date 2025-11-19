classdef ClusterConjugatePrior_t < ClusterConjugatePrior
    %CLUSTERCONJUGATEPRIOR Summary of this class goes here
    %   Detailed explanation goes here
    properties
        nu
        weights
        A
        sigma2
        
        ksai0
        tau0
%         CachedPosterior = []
    end
    
    methods
        function obj = ClusterConjugatePrior_t(ids, X, ts, params)
            %CLUSTER_DELTAMAP Construct an instance of this class
            %   Detailed explanation goes here
            obj@ClusterConjugatePrior(ids, X, ts, params);
            
            weights = parse_param(params, 'weights', []);
            obj.weights = weights;   
            obj.nu = parse_param(params, 'nu', []);
            
            obj.ksai0 = parse_param(params, 'ksai0', 1e-3);
            obj.tau0 = parse_param(params, 'tau0', 1e-3);            
            
            suffStats = calcSufficientStats_t(obj, obj.X, obj.V, weights);
            obj.SuffStats = suffStats;        
        end
        
        
        function addData(obj, id, x, t, w)
            % x and t and w are a cell array of multiple data
            obj.Ids = [obj.Ids, id];
            obj.X = [obj.X, x];
            obj.Ts = [obj.Ts, t];
            obj.weights = [obj.weights, w];
            
            [v, delta_t] = calculate_velocity(x, t);
            obj.V = [obj.V, v];
            obj.Delta_t = [obj.Delta_t, delta_t];
            
            % add sufficient statistics
            suffStats = obj.calcSufficientStats_t(x, v, w);
            obj.addSuffStats(suffStats);
        end
                
        function w = removeData(obj, id)
            % id should be a single data to remove
            k = find(obj.Ids == id);
            w = obj.weights(k);
            suffStats = obj.calcSufficientStats_t(obj.X(k), obj.V(k), obj.weights(k));
            
            selected = (obj.Ids ~= id);
            
            obj.Ids = obj.Ids(selected);
            obj.X = obj.X(selected);
            obj.Ts = obj.Ts(selected);
            obj.V = obj.V(selected);
            obj.Delta_t = obj.Delta_t(selected);
            obj.weights = obj.weights(selected);
            
            % remove sufficient stats
            obj.removeSuffStats(suffStats);
            
%             if length(obj.X) <= 1
%                 error('There is only one point left in the cluster');
%             end
        end
                
        function suffStats = calcSufficientStats_t(obj, X, V, weights)
            D = obj.Dim;
            XX1_new = zeros(D, D);
            VX1_new = zeros(D, D);
            V1V_new = 0;
            TiMinus1_new = 0;
            for i = 1:length(X)
                x1 = X{i}(:,1:end-1);
                x1w = x1 .* repmat(weights{i}, [D,1]);
                XX1_new = XX1_new + x1w * x1';
                VX1_new = VX1_new + V{i} * x1w';
                v1 = V{i} .* repmat(weights{i}, [D,1]);
                V1V_new = V1V_new + sum(sum(v1.*V{i}));
                TiMinus1_new = TiMinus1_new + size(V{i},2);
            end
            
            suffStats = [];
            suffStats.XX1 = XX1_new;
            suffStats.VX1 = VX1_new;
            suffStats.V1V = V1V_new;
            suffStats.TiMinus1 = TiMinus1_new;

        end
        
        function suffStats = calcSufficientStats(obj, X, V)
            suffStats = [];
        end

        
        function val = calcLogPredictive(obj, x, t, w)
            % calculate posterior parameters based on sufficient statistics
            [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParams(obj);
            
            % calculate current data point parameters
            [v, delta_t] = calculate_velocity(x, t);
            suffStats = calcSufficientStats_t(obj, x, v, w);
                        
            [d, Lambda, mu, eps] = calcLikelihoodParams(obj, suffStats);
            
            val = calcLogPredictiveByThm(obj, nu_p, ...
                Lambda_p, mu_p, kappa_p, Lambda, mu, eps, d, x, delta_t);
        end
        
        
        function [A, sigma2, A_mean, sigma2_mean] = sampleAsAndSigmas(obj)
            [A, sigma2, A_mean, sigma2_mean] = ...
                sampleAsAndSigmas@ClusterConjugatePrior(obj);
            obj.A = A;
            obj.sigma2 = sigma2;
        end
        
        function weights = sampleWeights(obj)
            D = obj.Dim;
            
            assert(length(obj.weights) == length(obj.X));
            
            for i = 1:length(obj.X)
                tmp = sum((obj.V{i} - obj.A * obj.X{i}(:,1:end-1)).^2, 1);
                rate = (tmp / obj.sigma2 + obj.nu) / 2;
                obj.weights{i} = gamrnd((D+obj.nu)/2, 1./rate);
            end
            
            suffStats = obj.calcSufficientStats_t(obj.X, obj.V, obj.weights);
            obj.SuffStats = suffStats;   
            
            weights = obj.weights;
        end
        
        function nu = sampleNu(obj)
            ksai1 = obj.SuffStats.TiMinus1 + obj.ksai0;
            tau1 = obj.tau0;
            for i = 1:length(obj.weights)
                tau1 = tau1 + sum((log(obj.weights{i}) - obj.weights{i}) / 2);
            end
            nu = sample_T_DOF(ksai1, tau1, 1);
            obj.nu = nu;
            
        end
    end

    methods(Static)
    end
end

