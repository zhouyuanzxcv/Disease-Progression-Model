function X = sample_dirichlet(a)
%SAMPLE_DIRICHLET 
% fastest method that generalize method used to sample from
% the BETA distribuition: Draw x1,...,xk from independent gamma 
% distribuitions with common scale and shape parameters a1,...,ak, 
% and for each j let rj=xj/(x1+...+xk).
% Input:
%   a - 1 x K alpha for dirichlet

K = size(a,2);
N = 1;

X = zeros(N,K);
for k = 1:K
    X(:,k) = gamrnd(a(k),1,N,1); % generates N random variables 
end                                   % with gamma distribution
X = X./repmat(sum(X,2),[1 K]);

% For more details see "Bayesian Data Analysis" Appendix A

end

