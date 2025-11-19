function E = calc_basis_for_centrosym_mat(D)
% construct basis matrix E for centrosymmetric matrix with size D x D
% Output:
%   E - basis matrix (size: D^2 x ceil(D^2/2)) such that vec(A) = E*a.
DOF = ceil((D^2)/2);
E_duplicate = zeros(D^2, D^2);
ind = 1;
for i = 1:D
    for j = 1:D
        tmp = zeros(D,D);
        tmp(i,j) = 1;
        tmp(D-i+1,D-j+1) = 1;
        E_duplicate(:,ind) = tmp(:);
        ind = ind + 1;
    end
end
E = unique(E_duplicate', 'rows')';
assert(size(E,2) == DOF);
end
