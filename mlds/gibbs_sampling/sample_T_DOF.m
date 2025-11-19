function nu = sample_T_DOF(ksai0, tau0, num_samples, domain)
%SAMPLE_T_DOF Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    domain = [1e-12, inf];
end

if nargin < 1
    test_sample_T_DOF
else
    nu = sample_T_DOF_impl(ksai0, tau0, num_samples, domain);
end

end

function nus = sample_T_DOF_impl(ksai0, tau0, num_samples, domain)
func = @logp1;
a = domain(1) + 1e-3;

% since the domain goes to infinity, logp(b) - logp(b-eps) should be less
% than 0
if isinf(domain(2))
    b = 100;
    for i = 1:10
        if func(b, ksai0, tau0, 0) - func(b - 1e-3, ksai0, tau0, 0) >= 0 
            b = b*10;
        else
            break;
        end
    end
else
    b = domain(2) - 1e-3;
end

% f_ab = func(linspace(a,b,10), ksai0, tau0, 0);

nus = ars(func, a, b, domain, num_samples, ksai0, tau0, 0);
end

function val = logp1(v, ksai0, tau0, const)
v1 = v/2;
val = ksai0 * (v1 .* log(v1) - gammaln(v1)) + v * tau0; 
val = val + const; % normalize the value to avoid too negative values
end

function test_sample_T_DOF
w = gamrnd(0.5, 1./(0.5*ones(1, 100)));
ksai0 = length(w);
tau0 = sum((log(w) - w)/2);
func = @logp1;
domain = [1e-12, inf];

samples = sample_T_DOF_impl(ksai0, tau0, 100, domain);


figure(3); clf;

subplot(1,2,1);
range = [1e-12 10];
x = (range(1):(range(2)-range(1))/(1000-1):range(2));
y = exp(feval(func,x,ksai0,tau0,0));
plot(x, y);
axis([range(1) range(2) 0 max(y)]);
title('f(x)');

subplot(1,2,2);

[N x] = hist(samples,100);
bar(x, N);
axis([ range(1) range(2) 0 max(N)]);
title('samples from f(x) (by ARS)');
end