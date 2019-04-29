function corrected = negative_correction(raw,sigma)

% cumulative normal from 0 to infinity of measurement -x, with std. dev. sigma
% expectation of one-sided truncation:
phi = @(x)(1/(sqrt(2*pi))*exp(-0.5*(x.^2)));
Phi = @(x)(0.5*(1+erf(x/sqrt(2))));

a = 0;
alpha = (a-raw)/sigma;
Z = 1-Phi(alpha);
corrected = raw + sigma*phi(alpha)./Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's try to compute this another way:
% Go up to k sigmas, see what the probability of it coming from this value is, assuming an even log distribution from 1-4sigma
bound = 4*sigma;
range = -bound:2*bound;
expectation = zeros(size(range));
for i=1:numel(range);
    expectation(i) = 0;
    cum_prob = 0;
    for j=0:0.01:log10(bound*10)
        prob = p_gaussian(range(i)-0.5,range(i)+0.5,10^(j+0.005),sigma);
        cum_prob = cum_prob + prob;
        expectation(i) = expectation(i) + (10.^(j+0.005))*prob;
    end
    expectation(i) = expectation(i)/cum_prob;
end

figure;
plot(range,expectation);

figure;
plot(range,log10(range)-log10(expectation))

rescales = (raw<max(range) & raw>min(range));
underflow = raw<=min(range);
corrected = raw;
corrected(underflow) = expectation(1);
corrected(rescales) = expectation(bound+ceil(0.5+raw(rescales)));


