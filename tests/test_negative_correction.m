% make a synthetic data set:
sigma = 100;
n_events = 1e5;
noise = randn(n_events,1)*sigma;
log_range = [0 5];
baseline = 10.^(log_range(1) + (rand(n_events,1)*(log_range(2)-log_range(1))));
cumulative = baseline+noise;

figure;
hist(noise,1e2);

figure;
hist(log10(abs(noise)),1e2);

figure;
hist(log10(baseline),1e2);

figure; hist(log10(cumulative(cumulative>0)),1e2);

figure; hist(cumulative(cumulative<800),1e2);


corrected = negative_correction(cumulative,sigma);
figure;
hist(log10(corrected),1e2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can we decide what an ideal corrective plot would be?
% First X% into first bin, etc...
