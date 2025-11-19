function X = fmin_centrosym(E,F)
%FMIN_CENTROSYM Summary of this function goes here
%   Detailed explanation goes here
D = size(E,1);

J = fliplr(eye(D));
invE = inv(E);
FinvE = F*invE;

tmp = (J*FinvE - FinvE*J)*inv(E*J*invE + J);

% calculate X
X = FinvE + tmp;

end

