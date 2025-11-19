function val = log_gamma(X, a, b)
val = a*log(b) - gammaln(a) + (a-1)*log(X) - b*X;
end


