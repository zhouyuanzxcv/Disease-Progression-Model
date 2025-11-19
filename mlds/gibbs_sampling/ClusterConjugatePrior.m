classdef ClusterConjugatePrior < TimeSeriesCluster
    %CLUSTERCONJUGATEPRIOR Summary of this class goes here
    %   Detailed explanation goes here
    properties
        % basis for the centrosymmetric matrices
        E
        
        % Parameters
        mu0
        Lambda0
        nu0
        kappa0
        
        use_centrosym
        
        
%         CachedPosterior = []
    end
    
    methods
        function obj = ClusterConjugatePrior(ids, X, ts, params)
            %CLUSTER_DELTAMAP Construct an instance of this class
            %   Detailed explanation goes here
            if ~isempty(X)
                D = size(X{1}, 1);
            else
                D = params.dim;
            end
            obj@TimeSeriesCluster(ids, X, ts, D, params);
            
            obj.use_centrosym = parse_param(params, 'use_centrosym', 1);
            if obj.use_centrosym
                obj.E = calc_basis_for_centrosym_mat(D);
            else
                E = parse_param(params, 'E_for_A', []);
                if isempty(E)
                    obj.E = eye(D^2);
                else
                    obj.E = E;
                end
            end
            DOF = size(obj.E, 2);
            
            obj.mu0 = parse_param(params, 'mu0', zeros(DOF,1));
            obj.Lambda0 = parse_param(params, 'Lambda0', eye(DOF)*1e-8);
            obj.nu0 = parse_param(params, 'nu0', 1e-3);
            obj.kappa0 = parse_param(params, 'kappa0', 1e-3);
            
        end
        
        function val = calcLogPredictive(obj, x, t)
            % calculate posterior parameters based on sufficient statistics
            [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParams(obj);
            
            % calculate current data point parameters
            [v, delta_t] = calculate_velocity(x, t);
            suffStats = calcSufficientStats(obj, x, v);
                        
            [d, Lambda, mu, eps] = calcLikelihoodParams(obj, suffStats);
            
            val = calcLogPredictiveByThm(obj, nu_p, ...
                Lambda_p, mu_p, kappa_p, Lambda, mu, eps, d, x, delta_t);
        end
        
        function [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParams(obj)
            [sum_d, sum_Lambda, sum_mu, sum_eps] = calcLikelihoodParams(...
                obj, obj.SuffStats);
            
            [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParamsByThm(obj, ...
                sum_d, sum_Lambda, sum_mu, sum_eps);
        end
        
        function [sum_d, sum_Lambda, sum_mu, sum_eps] = calcLikelihoodParams(...
                obj, suffStats)
            D = obj.Dim;

            sum_d = suffStats.TiMinus1 * D;
            sum_Lambda = obj.E' * kron(suffStats.XX1, eye(D)) * obj.E;
            sum_mu = obj.E' * suffStats.VX1(:);
            sum_eps = suffStats.V1V;
        end
        
        function [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParamsByThm(...
                obj, sum_d, sum_Lambda, sum_mu, sum_eps)
            % following Theorem 1
            nu_p = obj.nu0 + sum_d/2;
            Lambda_p = obj.Lambda0 + sum_Lambda;
            mu_p = Lambda_p \ (obj.Lambda0 * obj.mu0 + sum_mu);
            kappa_p = obj.kappa0 + 0.5 * (obj.mu0'*obj.Lambda0*obj.mu0 + ...
                sum_eps - mu_p'*Lambda_p*mu_p);
        end
        
        function [A, sigma2, A_mean, sigma2_mean] = sampleAsAndSigmas(obj)
            D = obj.Dim;
            [nu_p, Lambda_p, mu_p, kappa_p] = calcPosteriorParams(obj);
            
            % when there is only 1 sample, also throw error sometimes
            if ~isempty(obj.X)
                [a, sigma2] = sample_NIG(mu_p', Lambda_p, nu_p, kappa_p);
            else
                % given uninformative prior and an empty cluster
                a = nan(size(mu_p'));
                sigma2 = nan;
            end

            A = obj.E * a';
            A = reshape(A, [D D]);
            
            % set the expectation of NIG
            sigma2_mean = kappa_p / (nu_p - 1);
            
            A_mean = obj.E * mu_p;
            A_mean = reshape(A_mean, [D D]);
        end
        
        function val = calcLogPredictiveByThm(obj, nu_p, Lambda_p, mu_p, ...
                kappa_p, Lambda, mu, eps, d, x, delta_t)
            % following corollary
            
            % calculate Q
            tmp = Lambda_p*mu_p + mu;
            Q = mu_p'*Lambda_p*mu_p + eps - tmp' * ((Lambda_p + Lambda) \ tmp);
            
            % calculate log(p)
            % ignore log(h) since it is a common scalar in p(x_i|x_-i,k,
            % \beta) and p(x_i|\beta)
%             val = logmvn(x{1}(:,1)',obj.Params.mu_first,obj.Params.Sigma_first) ...
%                 - (d/2) * log(2*pi) - obj.Dim * sum(log(delta_t{1}));
            val = 0.5*logdet(Lambda_p) - 0.5*logdet(Lambda_p + Lambda) + ...
                gammaln(nu_p + d/2) - gammaln(nu_p) - (d/2)*log(kappa_p);
            val = val - (nu_p + d/2) * log(1 + (1/(2*kappa_p))*Q);
        end
    end

    methods(Static)
    end
end

