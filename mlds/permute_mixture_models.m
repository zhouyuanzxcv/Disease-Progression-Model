function models = permute_mixture_models(models)
%PERMUTE_MIXTURE_MODELS Summary of this function goes here
%   Detailed explanation goes here
[D,D,K] = size(models(1).As);
N = length(models);
ind = randperm(N,1);
template = reshape(models(ind).As, [D*D K])';
for i = 1:N
    model = models(i);
    As1 = reshape(model.As, [D*D K])';
    [P, As2] = permute_endmembers(template, As1);
    model.As = reshape(As2', [D D K]);
    model.sigmas = (P * model.sigmas')';
    model.pis = (P * model.pis')';
    model.mus_first = P * model.mus_first;
    model.Sigmas_first = reshape((P*reshape(model.Sigmas_first, [D*D K])')', [D D K]);
    models(i) = model;
end

end

