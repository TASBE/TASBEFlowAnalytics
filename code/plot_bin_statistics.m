% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_bin_statistics(sampleresults)

n_inductions = numel(sampleresults);
hues = (1:n_inductions)/n_inductions;

cfp_units = '';

stemName = TASBEConfig.get('OutputSettings.StemName');
directory = TASBEConfig.get('plots.plotPath');

%%% Bin count plots:
% Counts by CFP level:
maxcount = 1e1;
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
for i=1:n_inductions
    replicates = sampleresults{i};
    numReplicates = numel(replicates);
    for j=1:numReplicates,
        counts = replicates{j}.BinCounts;
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
            loglog(bin_centers(start:e),counts(start:e),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
            start = max(start+1, e+1);
        end
        % loglog(get_bin_centers(bins),counts,'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;

        %%% Gaussian Mixture model
        if(~isempty(replicates{j}.PlasmidModel))
            multiplier = sum(counts)*log10(bin_widths);
            fp_dist = get_fp_dist(replicates{j}.PlasmidModel);
            model = gmm_pdf(fp_dist, log10(bin_centers)')*multiplier;
            
            plot(bin_centers,model,'--','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
        %%% Old distribution model
    %         [fullmodel pbins] = CFP_distribution_model(replicates{j}.PlasmidModel);
    %         which = pbins>=min(get_bin_centers(bins)) & pbins<=max(get_bin_centers(bins));
    %         submodel = fullmodel(which)*sum(counts)/sum(fullmodel(which));
    %         loglog(pbins(which),submodel,'--','Color',hsv2rgb([hues(i) 1 0.9]));
        maxcount = max(maxcount,max(counts));
    end
end;
xlabel(['Constitutive ' clean_for_latex(cfp_units)]); ylabel('Count');
if(TASBEConfig.isSet('OutputSettings.FixedXAxis')), xlim(TASBEConfig.get('OutputSettings.FixedXAxis')); end;
if(TASBEConfig.isSet('OutputSettings.FixedYAxis')), ylim(TASBEConfig.get('OutputSettings.FixedYAxis')); else ylim([1e0 10.^(ceil(log10(maxcount)))]); end;
title([clean_for_latex(stemName),' bin counts, colored by inducer level']);
outputfig(h,[clean_for_latex(stemName),'-bincounts'],directory);

% Fraction active per bin:
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
for i=1:n_inductions
    replicates = sampleresults{i};
    numReplicates = numel(replicates);
    for j=1:numReplicates,
        analysisParam = replicates{j}.AnalysisParameters;
        bins = getBins(analysisParam);
        if ~isempty(replicates{j}.FractionActive),
            semilogx(get_bin_centers(bins),replicates{j}.FractionActive,'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
    end
end;
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel('Estimated Fraction Active');
ylim([-0.05 1.05]);
title([clean_for_latex(stemName),' estimated fraction of cells active, colored by inducer level']);
outputfig(h, [clean_for_latex(stemName),'-active'],directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasmid system is disabled, due to uncertainty about correctness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counts by estimated plasmid count:
% h = figure('PaperPosition',[1 1 5 3.66]);
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
