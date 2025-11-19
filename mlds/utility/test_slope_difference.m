function [p_diff] = test_slope_difference(xs,ys,coefs)
n1 = length(xs{1});
n2 = length(xs{2});

syx2_1 = calc_s_sqr(coefs{1}, xs{1}, ys{1});
syx2_2 = calc_s_sqr(coefs{2}, xs{2}, ys{2});

%% test if unequal residual variances
if syx2_1 > syx2_2
    ratio = syx2_1 / syx2_2;
    p = fcdf(ratio, n1-2, n2-2);
else
    ratio = syx2_2 / syx2_1;
    p = fcdf(ratio, n2-2, n1-2);
end

p = 2*(1-p);

Qx_1 = calc_Q(xs{1});
Qx_2 = calc_Q(xs{2});

if p < 0.05 % unequal residual variances
    z = abs(coefs{1}(1) - coefs{2}(1)) / sqrt(syx2_1 / Qx_1 + syx2_2 / Qx_2);
    p_diff = normcdf(-z)*2;
else % equal residual variances
    tmp1 = abs(coefs{1}(1) - coefs{2}(1));
    tmp2 = sqrt((syx2_1*(n1-2) + syx2_2*(n2-2)) / (n1+n2-4) * (1/Qx_1 + 1/Qx_2));
    t = tmp1 / tmp2;
    p_diff = tcdf(-t, n1+n2-4)*2;
end

end

% function p = calc_p_against_0(b, syx2, Qx, df)
% sbyx2 = syx2 / Qx;
% t = abs(b) / sqrt(sbyx2);
% p = tcdf(-t, df)*2;
% 
% end

function Q = calc_Q(x)
Q = sum((x - mean(x)).^2);
end

function s_sqr = calc_s_sqr(coef, x, y)
n = length(x);
y1h = x*coef(1) + coef(2);
s_sqr = sum((y - y1h).^2) / (n - 2);

end