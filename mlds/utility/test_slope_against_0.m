function p_against_0 = test_slope_against_0(x,y,coef)
% coef = [slope, intercept]
n = length(x);

syx2 = calc_s_sqr(coef, x, y);

Qx = calc_Q(x);

p_against_0 = calc_p_against_0(coef(1), syx2, Qx, n - 2);

end

function p = calc_p_against_0(b, syx2, Qx, df)
sbyx2 = syx2 / Qx;
t = abs(b) / sqrt(sbyx2);
p = tcdf(-t, df)*2;

end

function Q = calc_Q(x)
Q = sum((x - mean(x)).^2);
end

function s_sqr = calc_s_sqr(coef, x, y)
n = length(x);
y1h = x*coef(1) + coef(2);
s_sqr = sum((y - y1h).^2) / (n - 2);

end
