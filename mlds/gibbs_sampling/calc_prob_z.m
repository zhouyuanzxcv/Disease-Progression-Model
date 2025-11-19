function r_ik = calc_prob_z(zs)
[num_iter,N] = size(zs);
K = max(zs(:));
if K == 1 % only 1 cluster
    r_ik = ones(N,1);
    return;
end

zs1 = zeros(num_iter,N,K);

for k = 1:K
    zs1(:,:,k) = double(zs == k);
end

r_ik = squeeze(mean(zs1,1));

end


