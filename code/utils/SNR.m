function SNR_value = SNR(geomean1, geomean2, geostd1, geostd2)

% TODO: need to handle sigma better --- this is an approximation for when
% the stds are quite similar to one another
sigma = zeros(size(geostd1));
sigma(:) = geomean([geostd1(:) geostd2(:)]')';
SNR_value = 20 * log10( abs(log10(geomean1./geomean2))./(2*log10(sigma)));
