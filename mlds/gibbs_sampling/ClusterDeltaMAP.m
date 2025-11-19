classdef ClusterDeltaMAP < TimeSeriesCluster
    %CLUSTER_DELTAMAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DeltaMAP_InternalState = []
        
        % sigma0 is a parameter
        Params 
    end
    
    methods
        function obj = ClusterDeltaMAP(ids, X, ts, params)
            %CLUSTER_DELTAMAP Construct an instance of this class
            %   Detailed explanation goes here
            obj@TimeSeriesCluster(ids, X, ts, size(X{1},1));
            obj.Params = params;
            
            sigma0 = parse_param(params, 'sigma0', 1e4);
            obj.Params.sigma0_sqr = sigma0^2;
        end
        
        function val = calcLogPredictive(obj, x, t)
            [A, var] = estimateMAP(obj);
            val = calc_logprob_X_given_theta(x, t, A, sqrt(var), ...
                obj.Params.mu_first, obj.Params.Sigma_first);
        end
        
        function [A, sigma2] = sampleAsAndSigmas(obj)
            [A, sigma2] = estimateMAP(obj);
        end
                
        function [A, var] = estimateMAP(obj)
            D = obj.Dim;
            if ~isempty(obj.DeltaMAP_InternalState)
                %reuse last estimator if current sufficient stats are unchanged
                oldXX1 = obj.DeltaMAP_InternalState.oldXX1;
                oldVX1 = obj.DeltaMAP_InternalState.oldVX1;
                oldV1V = obj.DeltaMAP_InternalState.oldV1V;
                oldA = obj.DeltaMAP_InternalState.oldA;
                oldVar = obj.DeltaMAP_InternalState.oldVar;
            else
                % otherwise, initialize A and sigma 
                oldXX1 = nan(D, D);
                oldVX1 = nan(D, D);
                oldV1V = nan;
                oldA = fmin_centrosym(obj.SuffStats.XX1, obj.SuffStats.VX1);
                oldVar = obj.estimateSigma(oldA);
            end
            
            if all(all(oldXX1 == obj.SuffStats.XX1)) && ...
                    all(all(oldVX1 == obj.SuffStats.VX1)) && ...
                    oldV1V == obj.SuffStats.V1V
                % if the sufficient statistics is not changed, use the last
                % estimated A and variance
                A = obj.DeltaMAP_InternalState.oldA;
                var = obj.DeltaMAP_InternalState.oldVar;
            else
                % otherwise, estimate A and variance by block coordinate
                % descent
                [A, var] = obj.estimateMAP_bcd(oldA, oldVar);
                
                obj.DeltaMAP_InternalState.oldXX1 = obj.SuffStats.XX1;
                obj.DeltaMAP_InternalState.oldVX1 = obj.SuffStats.VX1;
                obj.DeltaMAP_InternalState.oldV1V = obj.SuffStats.V1V;
                obj.DeltaMAP_InternalState.oldA = A;
                obj.DeltaMAP_InternalState.oldVar = var;
            end
        end
        
        function [A, var] = estimateMAP_bcd(obj, A_init, var_init)
            A = A_init;
            var = var_init;
            D = size(A,1);
            maxA = max(abs(A(:)));
            
            oldA = A;
            for iter = 1:20
                E = obj.SuffStats.XX1 / var + eye(D) / obj.Params.sigma0_sqr;
                F = obj.SuffStats.VX1 / var;
                A = fmin_centrosym(E, F);
                var = estimateSigma(obj, A);
                
                if max(abs(A(:) - oldA(:))) < maxA * 1e-4
                    break;
                end
                
                oldA = A;
            end
        end
        
        function var = estimateSigma(obj, A)
            D = size(A,1);
            var = obj.SuffStats.V1V - 2*obj.SuffStats.VX1(:)'*A(:) + ...
                    A(:)'*kron(obj.SuffStats.XX1, eye(D))*A(:);
            var = var / (D * obj.SuffStats.TiMinus1 + 2);
        end

    end
end

