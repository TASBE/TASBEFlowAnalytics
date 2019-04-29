% PLOT_BATCH_HISTOGRAMS creates the bincount plots for batch analysis
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_batch_histograms(results,sampleresults,CM,linespecs)
% Elements of linespecs can either be LineSpecs, or ColorSpecs (e.g. three-element
% 0..1 vectors representing RGB color values); currently, only single-letter 
% color linespecs are properly handled.
figsize = TASBEConfig.get('OutputSettings.FigureSize');


% Obtain channel names in order to generate linespecs
channels = getChannelNames(sampleresults{1}{1}.AnalysisParameters); % channel names are same across conditions and replicates
if (~exist('linespecs', 'var'))
    % Build linespecs from sampleresults
    linespecs = cell(numel(channels),1);
    for i=1:numel(channels)
        if isempty(getLineSpec(channel_named(CM, channels{i})))
            linespecs{i} = 'k';
            TASBESession.warn('plot_batch_histograms','NoLineSpecs','Linespec for channel %s not found. Defaulting to black.', channels{i});
        else
            linespecs{i} = getLineSpec(channel_named(CM, channels{i}));
        end
    end
end
% Build unit names for X axis label
unitnames = {};
unitlegend = [];
for i=1:numel(channels)
    units = getUnits(channel_named(CM, channels{i}));
    if isempty(find(cellfun(@(u)(strcmp(u,units)),unitnames),1)) % if name is new
        unitnames{end+1} = units;
        if ~isempty(unitlegend), unitlegend = [unitlegend ' / ']; end;
        unitlegend = [unitlegend clean_for_latex(units)];
    end
end


if numel(linespecs) ~= numel(getChannelNames(sampleresults{1}{1}.AnalysisParameters))
    TASBESession.error('plot_batch_histograms','LineSpecDimensionMismatch', 'Size of linespecs does not match with number of channels');
end

n_conditions = size(sampleresults,1);
n_colors = numel(linespecs);

fprintf('Plotting histograms');

% Create legendentries
legendentries = getChannelNames(sampleresults{1}{1}.AnalysisParameters);

% one bincount plot per condition
maxcount = 1e1;
for i=1:n_conditions
    lines = [];
    % TODO: this should really be using a standard number
    h = figure('PaperPosition',[1 1 figsize]);
    set(h,'visible','off');
    bin_centers = results{i}.bincenters;
    for k=1:n_colors
        replicates = sampleresults{i};
        numReplicates = numel(replicates);
        for j=1:numReplicates
            counts = replicates{j}.BinCounts;
            ls = linespecs{k};
            isolates = isolated_points(counts(:,k),1);
            if(ischar(ls) && length(ls)==1 && length(findstr(ls, 'rgbcmykw')) == 1)
                line1 = loglog(bin_centers,counts(:,k),ls); hold on;
                loglog(bin_centers(isolates),counts(isolates,k),['+' ls]); % add isolated points with markers
            else
                line1 = loglog(bin_centers,counts(:,k),'Color', ls); hold on;
                loglog(bin_centers(isolates),counts(isolates,k),'+','Color',ls); % add isolated points with markers
            end
            
            if j == 1
                lines(end+1) = line1;
            end
        end
        maxcount = max(maxcount,max(max(counts)));
    end
    
    
    for j=1:numReplicates
        for k=1:n_colors
            ls = linespecs{k};
            if(ischar(ls) && length(ls)==1 && length(findstr(ls, 'rgbcmykw')) == 1)
                loglog([results{i}.means(k) results{i}.means(k)],[1 maxcount],[ls '--']); hold on;
            else
                loglog([results{i}.means(k) results{i}.means(k)],[1 maxcount], 'Color', ls, 'LineStyle', '--'); hold on;
            end
        end
    end
    
    xlabel(unitlegend); ylabel('Count');
    legend(lines, legendentries,'Location','Best');
    if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
    if(TASBEConfig.isSet('OutputSettings.FixedHistogramAxis')), ylim(TASBEConfig.get('OutputSettings.FixedHistogramAxis')); else ylim([1e0 10.^(ceil(log10(maxcount)))]); end

    title([TASBEConfig.get('OutputSettings.StemName') ' ' clean_for_latex(results{i}.condition) ' bin counts, by color']);
    
    if(strcmp(TASBEConfig.get('OutputSettings.StemName'), ' ') || strcmp(TASBEConfig.get('OutputSettings.StemName'), ''))
        outputfig(h,[results{i}.condition '-bincounts'],TASBEConfig.get('plots.plotPath'));
    else
        outputfig(h,[TASBEConfig.get('OutputSettings.StemName') '-' results{i}.condition '-bincounts'],TASBEConfig.get('plots.plotPath'));
    end
    
    fprintf('.');
end
fprintf('\n');
