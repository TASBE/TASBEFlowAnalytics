stem0312 = 'TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

% beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']); 
% data1 = fcs_scatter(beadfile,'FITC-A','Pacific Blue-A',0,[0 0; 6 6],1); % pure scatter
% data1 = fcs_scatter(beadfile,'FITC-A','Pacific Blue-A',0,[0 0; 6 6],1, false, [1 1]); % pure scatter, with linear set to true, doesn't seem to work
% data2 = fcs_scatter(beadfile,'FITC-A','Pacific Blue-A',1,[0 0; 6 6],1); % smoothed density plot

[rawdata, hdr, data] = fca_readfcs([stem0312 'Beads_P3.fcs']);

% let's pick one channel to deconvolve
channel = 'FITC-A'; % aka EYFP
sigma = 17.073972;
channel_data = getChannelData(hdr, data, channel);
[bins, counts] = getArithmeticBinning(channel_data, []);
counts_deconv = deconvolveData(counts, sigma);
% plot original counts with restored counts
plot(get_bin_centers(bins),counts,'b-'); hold on;
plot(get_bin_centers(bins),counts_deconv,'r-'); hold off;
xlabel('FITC a.u.'); ylabel('Count');
legend('counts with noise','deconvoluted counts')

% % create the scatter with x being FITC-A (deconvolved) and y being Pacific Blue-A
% smoothing = 10;
% % TODO: go from counts back to channel data
% % Error since counts_deconv not all integers, not sure what to do
% % channel_data_deconv = cell2mat(arrayfun(@(a, f) repmat(a, [1 f]), get_bin_centers(bins), counts_deconv, 'UniformOutput', false));
% channel_data2 = getChannelData(hdr, data, 'Pacific Blue-A');
% [bins2, counts2] = getArithmeticBinning(channel_data2, bins);
% num_bins = numel(get_bin_centers(bins));
% % getLogBinning(channel_data, []);
% smoothhist2D([channel_data channel_data2], smoothing, [num_bins num_bins]);

function channel_data = getChannelData(hdr, data, channel_name)
    channel_data = [];
    % get part of data that corresponds to selected channel
    idx = 1;
    for name = {hdr.par(:).name}
        if(strcmp(name{1}, channel_name))
            channel_data = data(:, idx);
            display(size(channel_data)); % size of channel_data
            display(sum(channel_data<=0)); % num entries that are negative
        end
        idx = idx + 1;
    end
end

function [bins, counts] = getArithmeticBinning(channel_data, bins)
    if numel(bins) == 0
        % create arithmetic binning of channel_data
        bins = BinSequence(min(channel_data),1,max(channel_data),'arithmetic');
        counts = histcounts(channel_data,get_bin_edges(bins));
    else
        counts = histcounts(channel_data,get_bin_edges(bins));
    end
    % plot(get_bin_centers(bins), counts,'b-');
    % xlabel('FITC a.u.'); ylabel('Count');

    % visualize counts as an image
    % original_image = imagesc(counts);
    % title('Counts with noise');
end

function [bins, counts] = getLogBinning(channel_data, bins)
    if numel(bins) == 0
        % create log binning of channel_data
        rmin = min(channel_data);
        rmax = max(channel_data);
        bins = BinSequence((rmin-0.5),0.05,(rmax+0.5),'log_bins');
        counts = histcounts(channel_data,get_bin_edges(bins));
    else
        counts = histcounts(channel_data,get_bin_edges(bins));
    end
    % semilogx(get_bin_centers(bins), counts, 'b-');
    % xlim([1e0 1e5]);
    % xlabel('log10 FITC a.u.'); ylabel('Count');
end

function counts_deconv = deconvolveData(counts, sigma)
    % create Gaussian kernel with standard deviation from autofluorescence
    % model 
    PSF = fspecial('gaussian',[1 numel(counts)],sigma);
    % display(sum(PSF));
    % plot(PSF);
    % xlabel('Pixel number'); ylabel('Gaussian function output');

    % apply deconv method Richardson-Lucy (RL) deconvolution
    counts_deconv = deconvlucy(counts,PSF); % default number of iterations is 10
    % imshow(luc1);
    % title('Restored Image');
end
