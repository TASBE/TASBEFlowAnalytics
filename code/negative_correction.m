function [corrected,cum_prob] = negative_correction(raw,sigma)

% Building the model is expensive, takes time proportional to sigma, and can be done with an autofluorescence model
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's try to compute this another way:
% Go up to k sigmas, see what the probability of it coming from this value is, assuming an even log distribution from 1-4sigma
bound = 4*sigma;
range = -bound:8*bound;
inverse_grain = 0.01;
inverse_range = 0:inverse_grain:log10(bound*8);
cum_prob = zeros(numel(range),numel(inverse_range));
for i=1:numel(range);
    range_prob = zeros(size(inverse_range));
    for j=1:numel(inverse_range)
        range_prob(j) = p_gaussian(range(i)-0.5,range(i)+0.5,10^(inverse_range(j)-inverse_grain/2),sigma);
    end
    cum_prob(i,:) = cumsum(range_prob)/sum(range_prob); % normalized distribution
end
toc

% Applying is proportional to both sigma and number of events, and is fairly fast
tic

rescales = (raw<max(range) & raw>min(range));
underflow = raw<=min(range);
corrected = raw;
corrected(underflow) = 1; % 10^0

n_rescales = sum(rescales);
rescale_projection = rand(n_rescales,1);
old_rescale_values = raw(rescales);
rescale_probability = 1-((old_rescale_values/(max(range)/2))-1);
rescale_feathering = rand(n_rescales,1);
rescale_indices = bound+ceil(0.5+old_rescale_values);
rescale_values = nan(n_rescales,1);
rescale_dequantize = (rand(n_rescales,1)-1)*inverse_grain; % randomly select within the range to the mark
for i=1:n_rescales
    if rescale_feathering(i) < rescale_probability(i)
        rescale_values(i) = 10.^(rescale_dequantize(i) + inverse_range(find(rescale_projection(i)<cum_prob(rescale_indices(i),:),1)));
    else
        rescale_values(i) = old_rescale_values(i);
    end
end

corrected(rescales) = rescale_values;
toc

