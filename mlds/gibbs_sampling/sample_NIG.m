function [a, sigma2] = sample_NIG(mu, Lambda, nu, kappa)
%SAMPLE_NIG Generate random samples from normal inverse gamma
% Input:
%   mu - mean of normal (1 x D)
%   Lambda - precision matrix for normal (D x D)
%   nu - shape of inverse gamma
%   kappa - rate of inverse gamma
% Output:
%   a - x from the normal
%   sigma2 - squared sigma from inverse gamma
sigma2 = 1 / gamrnd(nu, 1/kappa);

cov1 = sigma2*inv(Lambda);

if any(isinf(cov1(:)))
    a = nan(size(mu));
else
    % check symmetry
    if max(max(abs(cov1 - cov1'))) < max(abs(cov1(:))) * 1e-6
        cov1 = (cov1 + cov1')/2;
    else
        error('covariance is not symmetric in sample_NIG.');
    end

    a = mvnrnd(mu, cov1);
end

end

