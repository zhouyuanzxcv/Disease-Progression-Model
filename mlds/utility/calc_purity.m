function purity = calc_purity(z1, z2)
%CALC_PURITY Summary of this function goes here
%   Detailed explanation goes here
K1 = max(z1);
K2 = max(z2);
N = length(z1);

for k = 1:K1
    inds1{k} = find(z1 == k);
end

for k = 1:K2
    inds2{k} = find(z2 == k);
end

max_szs = [];
for k = 1:K1
    inter_szs = [];
    for k2 = 1:K2
        inter_szs(k2) = length(intersect(inds1{k}, inds2{k2}));
    end
    [m,ind] = max(inter_szs);
    max_szs(k) = m;
end

purity = sum(max_szs) / N;

end

