% PLOT_BIN_STATISTICS_ALL creates bincount plots by CFP level and fraction
% active per bin for plusminus and transfer curve analysis with results from 
% all condition levels combined into one plot.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_bin_statistics_all(allsampleresults, batch_names, inducer_levels, index, channelName)
cfp_units = '';

stemName = TASBEConfig.get('OutputSettings.StemName');
directory = TASBEConfig.get('plots.plotPath');
figsize = TASBEConfig.get('OutputSettings.FigureSize');

% Setting up additional plot features 
n_var = numel(allsampleresults{1});
n_batch = numel(batch_names);
% Determine which dimension is varied by which line spec property (marker
% or color)
ref_marker_types = {'o', '*', '.', 'x'};
by_n_var = 0; % state variable
if n_var > n_batch
    hues = (1:n_var)./n_var;
    marker_types = ref_marker_types(1:n_batch);
    by_n_var = 1;
else
    marker_types = ref_marker_types(1:n_var);
    hues = (1:n_batch)./(n_batch);
end

% Create legendentries
legendentries = cell(0);
for j=1:n_batch
    if j == 1
        for k=1:n_var
            entrystr = num2str(inducer_levels(k));
            legendentries{end+1} = [char(batch_names{1}) ' ' entrystr]; % not pre-allocated because we don't know how many are valid
        end
    else
        entrystr = char(batch_names{j});
        legendentries{end+1} = entrystr;    
    end
end

%%% Bin count plots:
% Counts by CFP level:
maxcount = 1e1;
h = figure('PaperPosition',[1 1 figsize]);
set(h,'visible','off');
lines = [];
legendentries2 = legendentries;
legendentries2{end+1} = 'Gaussian Mixture model';
for k=1:n_batch
    sampleresults = allsampleresults{k};
    for i=1:n_var
        replicates = sampleresults{i};
        numReplicates = numel(replicates);
        for j=1:numReplicates
            counts = replicates{j}.Histograms(:,index);
            % counts = replicates{j}.BinCounts;
            analysisParam = replicates{j}.AnalysisParameters;
            bins = getBins(analysisParam);
            bin_centers = get_bin_centers(bins);
            bin_widths = get_bin_widths(bins);
            rep_units = getChannelUnits(analysisParam,'constitutive');
            if strcmp(cfp_units,''), cfp_units = rep_units;
            else if ~strcmp(cfp_units,rep_units), cfp_units = 'a.u.';
                end
            end

            start = 1;
            bin_size = numel(bin_centers);

            while (start <= bin_size)
                nanLoc = find(isnan(bin_centers(start:bin_size)), 1);
                if isempty(nanLoc)
                    e = bin_size;
                else
                    e = nanLoc - 1 + start - 1;
                end
                if by_n_var
                    line1 = loglog(bin_centers(start:e),counts(start:e),['-' marker_types{k}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
                    if(start == e), loglog(bin_centers(start),counts(start),marker_types{k},'Color',hsv2rgb([hues(i) 1 0.9])); hold on; end % make sure isolated points show
                    start = max(start+1, e+1);
                else
                    line1 = loglog(bin_centers(start:e),counts(start:e),['-' marker_types{i}],'Color',hsv2rgb([hues(k) 1 0.9])); hold on;
                    if(start == e), loglog(bin_centers(start),counts(start),marker_types{i},'Color',hsv2rgb([hues(k) 1 0.9])); hold on; end % make sure isolated points show
                    start = max(start+1, e+1);
                end
            end
            % loglog(get_bin_centers(bins),counts,'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;

            %%% Gaussian Mixture model
            if(~isempty(replicates{j}.PlasmidModel) && ~isempty(get_fp_dist(replicates{j}.PlasmidModel)))
                multiplier = sum(counts)*log10(bin_widths);
                fp_dist = get_fp_dist(replicates{j}.PlasmidModel);
                model = gmm_pdf(fp_dist, log10(bin_centers)')*multiplier;
                if by_n_var
                    line2 = plot(bin_centers,model,['--' marker_types{k}],'Color',hsv2rgb([hues(i) 1 0.9]),'MarkerSize', 3); hold on;
                else
                    line2 = plot(bin_centers,model,['--' marker_types{i}],'Color',hsv2rgb([hues(k) 1 0.9]),'MarkerSize', 3); hold on;
                end
            end
            %%% Old distribution model
        %         [fullmodel pbins] = CFP_distribution_model(replicates{j}.PlasmidModel);
        %         which = pbins>=min(get_bin_centers(bins)) & pbins<=max(get_bin_centers(bins));
        %         submodel = fullmodel(which)*sum(counts)/sum(fullmodel(which));
        %         loglog(pbins(which),submodel,'--','Color',hsv2rgb([hues(i) 1 0.9]));
            maxcount = max(maxcount,max(counts));

            % Determines which type of lines should be in the legend
            if j == 1 && (i == 1 || k == 1 )
                lines(end+1) = line1;
            end

            if i == n_var && j == numReplicates && k == n_batch
                lines(end+1) = line2;
            end
        end
    end
end
xlabel([channelName ' ' clean_for_latex(cfp_units)]); ylabel('Count');
legend(lines, legendentries2,'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedHistogramAxis')), ylim(TASBEConfig.get('OutputSettings.FixedHistogramAxis')); else ylim([1e0 10.^(ceil(log10(maxcount)))]); end
title(['All ' clean_for_latex(stemName),' bin counts, colored by inducer level']);
if(strcmp(clean_for_latex(stemName), ' ') || strcmp(clean_for_latex(stemName), ''))
    outputfig(h,['all-' channelName 'bincounts'],directory);
else
    outputfig(h,[clean_for_latex(stemName),'-all-' channelName '-bincounts'],directory);
end

% Fraction active per bin:
h = figure('PaperPosition',[1 1 figsize]);
set(h,'visible','off');
lines = [];
for k=1:n_batch
    sampleresults = allsampleresults{k};
    for i=1:n_var
        replicates = sampleresults{i};
        numReplicates = numel(replicates);
        for j=1:numReplicates
            analysisParam = replicates{j}.AnalysisParameters;
            bins = getBins(analysisParam);
            if ~isempty(replicates{j}.FractionActive)
                if by_n_var
                    line1 = semilogx(get_bin_centers(bins),replicates{j}.FractionActive,['-' marker_types{k}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
                else
                    line1 = semilogx(get_bin_centers(bins),replicates{j}.FractionActive,['-' marker_types{i}],'Color',hsv2rgb([hues(k) 1 0.9])); hold on;
                end
            end

            % Determines which type of lines should be in the legend
            if j == 1 && (i == 1 || k == 1 )
                lines(end+1) = line1;
            end
        end
    end
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel('Estimated Fraction Active');
legend(lines, legendentries,'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
ylim([-0.05 1.05]);
title(['All ' clean_for_latex(stemName),' estimated fraction of cells active, colored by inducer level']);
if(strcmp(clean_for_latex(stemName), ' ') || strcmp(clean_for_latex(stemName), ''))
    outputfig(h, 'all-active',directory);
else
    outputfig(h, [clean_for_latex(stemName),'-all-active'],directory);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasmid system is disabled, due to uncertainty about correctness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counts by estimated plasmid count:
% h = figure('PaperPosition',[1 1 figsize]);
% set(h,'visible','off');
% for i=1:n_inductions
%     replicates = sampleresults{i};
%     numReplicates = numel(replicates);
%     for j=1:numReplicates,
%         pe = replicates{j}.PlasmidEstimates;
%         counts = replicates{j}.BinCounts;
%         active = replicates{j}.FractionActive;
%         pe(active<0.9) = NaN;
%         
%         start = 1;
%         pe_size = size(pe, 1);
%         while (start <= pe_size)
%             nanLoc = find(isnan(pe(start:pe_size)), 1);
%             if isempty(nanLoc)
%                 e = pe_size;
%             else
%                 e = nanLoc - 1 + start - 1;
%             end
%             loglog(pe(start:e),counts(start:e),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
%             start = max(start+1, e+1);
%         end
%         
%         % loglog(replicates{j}.PlasmidEstimates,replicates{j}.BinCounts,'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
%     end
% end;
% set(gca,'XScale','log'); set(gca,'YScale','log');
% xlabel('Estimated Plasmid Count'); ylabel('Count');
% ylim([1e0 10.^(ceil(log10(maxcount)))]);
% title([TASBEConfig.get('OutputSettings.StemName'),' bin counts, colored by inducer level']);
% outputfig(h,[TASBEConfig.get('OutputSettings.StemName'),'-plasmid-bincounts'],TASBEConfig.get('OutputSettings.Directory'));
