function corrected = negative_correction(raw,sigma)

% cumulative normal from 0 to infinity of measurement -x, with std. dev. sigma
% expectation of one-sided truncation:
phi = @(x)(1/(sqrt(2*pi))*exp(-0.5*(x.^2)));
Phi = @(x)(0.5*(1+erf(x/sqrt(2))));

a = 0;
alpha = (a-raw)/sigma;
Z = 1-Phi(alpha);
corrected = raw + sigma*phi(alpha)./Z;
