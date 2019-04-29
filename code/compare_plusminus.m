% COMPARE_PLUSMINUS(p_results,m_results,no_input)
%    p_results and m_results are the ExperimentResults from the positive and negative conditions
%    if 'no_input' (optional) is set, then the input value will be set uniformly to the induction level
%
%    imeans, omeans are [bins x inductions x 2] arrays of input, output means
%    bincounts is a [bins x inductions x 2] array of bin counts for +, -
%    ratios is a [bins x inductions] array of +/- ratio
%    ostds, ostdofmeans are variances and variances of omeans
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function pm_results = compare_plusminus(p_results,m_results,no_input)
if(nargin<3), no_input = 0; end; % default to input, not only inducer

% Ensure results are comparable
p_AP = getAnalysisParameters(p_results);
m_AP = getAnalysisParameters(m_results);
p_bins = getBins(p_AP);
m_bins = getBins(m_AP);
p_inductions = getInducerLevelsToFiles(getExperiment(p_results),1);
m_inductions = getInducerLevelsToFiles(getExperiment(m_results),1);

if ((get_n_bins(p_bins) ~= get_n_bins(m_bins)) || ...
    numel(find(get_bin_edges(p_bins) ~= get_bin_edges(m_bins))) || ...
    numel(find(p_inductions ~= m_inductions)))
    TASBESession.error('TASBE:Analysis','PlusMinusMismatch','Plus/Minus experiment cannot compare datasets with different inductions or bins');
end

n_bins = get_n_bins(p_bins);
inductions = p_inductions;


% Obtain raw results
[p_omeans p_ostds p_ostdofmeans] = get_channel_results(p_results,'output');
if no_input==0, [p_imeans  p_istds p_istdofmeans]= get_channel_results(p_results,'input');
else p_imeans = ones(n_bins,1)*inductions; p_istds = zeros(n_bins,1); p_istdofmeans = zeros(n_bins,1);
end
p_active = getFractionActive(p_results);
p_bincounts = getBinCounts(p_results);
p_valid = p_active >= getMinFractionActive(p_AP) & ...
    p_bincounts >= getMinValidCount(p_AP);

[m_omeans m_ostds m_ostdofmeans] = get_channel_results(m_results,'output');
if no_input==0, [m_imeans m_istds m_istdofmeans] = get_channel_results(m_results,'input');
else m_imeans = ones(n_bins,1)*inductions; m_istds = zeros(n_bins,1); m_istdofmeans = zeros(n_bins,1);
end
m_active = getFractionActive(m_results);
m_bincounts = getBinCounts(m_results);
m_valid = m_active >= getMinFractionActive(m_AP) & ...
    m_bincounts >= getMinValidCount(m_AP);

% Zip together into output
imeans = zeros(n_bins,numel(inductions),2); 
omeans = imeans; ostds = omeans; ostdofmeans = omeans; bincounts = omeans; valid = omeans;
imeans(:,:,1) = p_imeans; imeans(:,:,2) = m_imeans;
istds(:,:,1) = p_istds; istds(:,:,2) = m_istds;
istdofmeans(:,:,1) = p_istdofmeans; istdofmeans(:,:,2) = m_istdofmeans;
omeans(:,:,1) = p_omeans; omeans(:,:,2) = m_omeans;
ostds(:,:,1) = p_ostds; ostds(:,:,2) = m_ostds;
ostdofmeans(:,:,1) = p_ostdofmeans; ostdofmeans(:,:,2) = m_ostdofmeans;
bincounts(:,:,1) = p_bincounts; bincounts(:,:,2) = m_bincounts;
valid(:,:,1) = p_valid; valid(:,:,2) = m_valid;

bothvalid = p_valid & m_valid;
pm_ratios = omeans(:,:,1)./omeans(:,:,2);
iSNR = SNR(imeans(:,:,1), imeans(:,:,2), istds(:,:,1), istds(:,:,2));
oSNR = SNR(omeans(:,:,1), omeans(:,:,2), ostds(:,:,1), ostds(:,:,2));

mean_ratio = zeros(size(pm_ratios,2),1);
for j=1:size(pm_ratios,2),
    mean_ratio(j) = mean(pm_ratios(bothvalid(:,j),j)); % get mean ratio for each condition
end


pm_results = PlusMinusResults(p_results, m_results, mean_ratio, bincounts, valid,...
    pm_ratios, imeans, istds, istdofmeans, omeans, ostds, ostdofmeans);

pm_results.InputSNR = iSNR;
pm_results.OutputSNR = oSNR;
