% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_plusminus_comparison(all_pm_results, batch_names)
% Obtaining TASBEConfig preferences for various graphs
step = TASBEConfig.get('OutputSettings.PlotEveryN');
ticks = TASBEConfig.get('OutputSettings.PlotTickMarks');
deviceName = TASBEConfig.get('OutputSettings.DeviceName');
directory = TASBEConfig.get('plots.plotPath');
stemName = TASBEConfig.get('OutputSettings.StemName');

% Setting up additional plot features 
variable = getInducerLevelsToFiles(getExperiment(all_pm_results{1}.PlusResults),1);
n_var = numel(variable);
n_comparisons = numel(all_pm_results);
% Determine which dimension is varied by which line spec property (marker
% or color)
ref_marker_types = {'o', '*', '.', 'x'};
by_n_var = 0; % state variable
if n_var > n_comparisons+1
    hues = (1:n_var)./n_var;
    marker_types = ref_marker_types(1:n_comparisons+1);
    by_n_var = 1;
else
    marker_types = ref_marker_types(1:n_var);
    hues = (1:n_comparisons+1)./(n_comparisons+1);
end

% Gather all of the AP, bin_centers, and units
all_AP = {};
all_bin_centers = {};
for i=1:n_comparisons
    AP = getAnalysisParameters(all_pm_results{i}.PlusResults);
    bin_centers = get_bin_centers(getBins(AP));
    all_AP{i} = AP;
    all_bin_centers{i} = bin_centers;
end
in_units = getChannelUnits(all_AP{1},'input');
out_units = getChannelUnits(all_AP{1},'output');
cfp_units = getChannelUnits(all_AP{1},'constitutive');

% Create legendentries
legendentries = cell(0);
for j=1:n_comparisons+1
    if j == 1
        for k=1:step:n_var
            which = all_pm_results{1}.Valid(:,k,1) & all_pm_results{1}.Valid(:,k,2);
            entrystr = num2str(variable(k));
            if numel(which)>0 && sum(which)>0
                legendentries{end+1} = [char(batch_names{1}) ' ' entrystr]; % not pre-allocated because we don't know how many are valid
            else
                warning('PlotPlusMinus:EmptyResults','No active results for %s: not enough data or bad component fit',entrystr);
            end
        end
    else
        entrystr = char(batch_names{j});
        legendentries{end+1} = entrystr;    
    end
end

% Legendentries for comparison graphs
comlegendentries = legendentries(1:end-1);

%%% I/O plots:
% Plain I/O plot:
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
add_lines = true;
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    % Plot the results for each batch_name and level
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
            line2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2),['-' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
            line2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2),['-' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9])); hold on;
        end
        
        % To avoid duplicating lines, only plot the last batch_name once
        if i ~= numel(all_pm_results)
            set(line2, 'visible', 'off');
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
            if i == numel(all_pm_results) && add_lines
                lines(end+1) = line2;
                add_lines = false;
            end
        end

    end
    
    % Plot error envelope if specified
    if TASBEConfig.get('plusminus.plotError')
        for j=1:step:n_var
            which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
            if by_n_var
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1).*pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2).*pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
            else
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1).*pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2).*pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
            end
            
            if i ~= numel(all_pm_results)
                set(line2_1, 'visible', 'off');
                set(line2_2, 'visible', 'off');
            end
        end
    end
    hold on;
end
xlabel(['IFP ' clean_for_latex(in_units)]); ylabel(['OFP ' clean_for_latex(out_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
legend(lines, legendentries,'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedOutputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedOutputAxis')); end
title(['Raw All ',clean_for_latex(stemName),' Transfer Curves']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-mean'],directory);

% normalized I/O plot
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
add_lines = true;
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    % Plot the results for each batch_name and level
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./all_bin_centers{i}(which)',['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
            line2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./all_bin_centers{i}(which)',['-' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./all_bin_centers{i}(which)',['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
            line2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./all_bin_centers{i}(which)',['-' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9])); hold on;
        end
        
        % To avoid duplicating lines, only plot the last batch_name once
        if i ~= numel(all_pm_results)
            set(line2, 'visible', 'off');
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
            if i == numel(all_pm_results) && add_lines
                lines(end+1) = line2;
                add_lines = false;
            end
        end

    end
    % Plot error envelope if specified
    if TASBEConfig.get('plusminus.plotError')
        for j=1:step:n_var
            which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
            if by_n_var
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./bin_centers(which)'.*pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./bin_centers(which)'./pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./bin_centers(which)'.*pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./bin_centers(which)'./pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
            else
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./bin_centers(which)'.*pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                loglog(pm_results.InMeans(which,j,1),pm_results.OutMeans(which,j,1)./bin_centers(which)'./pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./bin_centers(which)'.*pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(pm_results.InMeans(which,j,2),pm_results.OutMeans(which,j,2)./bin_centers(which)'./pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
            end
            
            if i ~= numel(all_pm_results)
                set(line2_1, 'visible', 'off');
                set(line2_2, 'visible', 'off');
            end
        end
    end
    hold on;
end
xlabel(['IFP ' clean_for_latex(in_units)]); ylabel(['OFP ' clean_for_latex(out_units) ' / CFP ' clean_for_latex(cfp_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
legend(lines, clean_for_latex(legendentries),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedNormalizedOutputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedNormalizedOutputAxis')); end
title(['All ' clean_for_latex(stemName),' Transfer Curves Normalized by CFP']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-mean-norm'],directory);

% IFP vs. CFP
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
add_lines = true;
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    % Plot the results for each batch_name and level
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
            line2 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2),['-' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
            line2 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2),['-' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9])); hold on;
        end
        
        % To avoid duplicating lines, only plot the last batch_name once
        if i ~= numel(all_pm_results)
            set(line2, 'visible', 'off');
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
            if i == numel(all_pm_results) && add_lines
                lines(end+1) = line2;
                add_lines = false;
            end
        end

    end
    
    % Plot error envelope if specified
    if TASBEConfig.get('plusminus.plotError')
        for j=1:step:n_var
            which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
            if by_n_var
                loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1).*pm_results.InStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1)./pm_results.InStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2).*pm_results.InStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2)./pm_results.InStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
            else
                loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1).*pm_results.InStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,1)./pm_results.InStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2).*pm_results.InStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(all_bin_centers{i}(which),pm_results.InMeans(which,j,2)./pm_results.InStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
            end
            
            if i ~= numel(all_pm_results)
                set(line2_1, 'visible', 'off');
                set(line2_2, 'visible', 'off');
            end
        end
    end
    hold on;
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel(['IFP ' clean_for_latex(out_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
legend(lines, clean_for_latex(legendentries),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end
title(['All ' clean_for_latex(stemName),' IFP vs. CFP']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-input-v-cfp'],directory);

% OFP vs. CFP
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
add_lines = true;
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    % Plot the results for each batch_name and level
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
            line2 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2),['-' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
            line2 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2),['-' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9])); hold on;
        end
        
        % To avoid duplicating lines, only plot the last batch_name once
        if i ~= numel(all_pm_results)
            set(line2, 'visible', 'off');
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
            if i == numel(all_pm_results) && add_lines
                lines(end+1) = line2;
                add_lines = false;
            end
        end

    end
    
    % Plot error envelope if specified
    if TASBEConfig.get('plusminus.plotError')
        for j=1:step:n_var
            which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
            if by_n_var
                loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1).*pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1)./pm_results.OutStandardDevs(which,j,1),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2).*pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2)./pm_results.OutStandardDevs(which,j,2),[':' marker_types{i+1}],'Color',hsv2rgb([hues(j) 1 0.9]), 'MarkerSize', 3);
            else
                loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1).*pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,1)./pm_results.OutStandardDevs(which,j,1),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9]), 'MarkerSize', 3);
                line2_1 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2).*pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
                line2_2 = loglog(all_bin_centers{i}(which),pm_results.OutMeans(which,j,2)./pm_results.OutStandardDevs(which,j,2),[':' marker_types{j}],'Color',hsv2rgb([hues(i+1) 1 0.9]), 'MarkerSize', 3);
            end
            
            if i ~= numel(all_pm_results)
                set(line2_1, 'visible', 'off');
                set(line2_2, 'visible', 'off');
            end
        end
    end
    hold on;
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel(['OFP ' clean_for_latex(out_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
legend(lines, clean_for_latex(legendentries),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedOutputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedOutputAxis')); end
title(['All ' clean_for_latex(stemName),' OFP vs. CFP']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-v-cfp'],directory);

% % Relative change in OFP vs. CFP
% Removed because it wasn't ever useful
% h = figure('PaperPosition',[1 1 5 3.66]);
% set(h,'visible','off');
% for i=1:step:n_var
%     which = find(pm_results.Valid(1:(end-1),i,1) & pm_results.Valid(1:(end-1),i,2) & ...
%         pm_results.Valid(2:end,i,1) & pm_results.Valid(2:end,i,2));
%     marginal_centers = (bin_centers(which) + bin_centers(which+1)) / 2;
%     p_ofp_difference = pm_results.OutMeans(which+1,i,1) ./ pm_results.OutMeans(which,i,1);
%     m_ofp_difference = pm_results.OutMeans(which+1,i,2) ./ pm_results.OutMeans(which,i,2);
%     cfp_difference = (bin_centers(which+1) ./ bin_centers(which));
%     semilogx(marginal_centers,p_ofp_difference./cfp_difference',[ptick '-'],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
% end
% for i=1:step:n_var
%     which = find(pm_results.Valid(1:(end-1),i,1) & pm_results.Valid(1:(end-1),i,2) & ...
%         pm_results.Valid(2:end,i,1) & pm_results.Valid(2:end,i,2));
%     marginal_centers = (bin_centers(which) + bin_centers(which+1)) / 2;
%     p_ofp_difference = pm_results.OutMeans(which+1,i,1) ./ pm_results.OutMeans(which,i,1);
%     m_ofp_difference = pm_results.OutMeans(which+1,i,2) ./ pm_results.OutMeans(which,i,2);
%     cfp_difference = (bin_centers(which+1) ./ bin_centers(which));
%     semilogx(marginal_centers,m_ofp_difference./cfp_difference',[ntick '--'],'Color',hsv2rgb([hues(i) 1 0.9]));
% end;
% xlabel(['CFP ' cfp_units]); ylabel(['OFP ' out_units ' / CFP ' cfp_units]);
% set(gca,'XScale','log'); set(gca,'YScale','linear');
% legend('Location','Best',pmlegendentries,'Minus');
% if(OutputSettings.FixedInputAxis), xlim(OutputSettings.FixedInputAxis); end;
% if(OutputSettings.FixedOutputAxis), ylim(OutputSettings.FixedOutputAxis); end;
% title([OutputSettings.StemName,' marginal change in OFP vs. CFP']);
% outputfig(h,[OutputSettings.StemName,'-',OutputSettings.DeviceName,'-marginal-ofp'],directory);
% 

% ratio plot
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = semilogx(all_bin_centers{i}(which),pm_results.Ratios(which,j),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = semilogx(all_bin_centers{i}(which),pm_results.Ratios(which,j),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
        end
    end
    hold on;
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel('Fold Activation');
set(gca,'XScale','log'); set(gca,'YScale','log');
legend(lines, clean_for_latex(comlegendentries),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedRatioAxis')), ylim(TASBEConfig.get('OutputSettings.FixedRatioAxis')); end
title(['All Ratios (in respect to ' batch_names{end} ') for ',clean_for_latex(stemName)]);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-ratios'],directory);

% SNR plots
legendentries2 = legendentries;
for i=1:numel(legendentries)
    legendentries2{i} = [legendentries{i} ' output SNR'];
end
legendentries2{end} = 'input SNR';
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = semilogx(all_bin_centers{i}(which),pm_results.OutputSNR(which,j),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = semilogx(all_bin_centers{i}(which),pm_results.OutputSNR(which,j),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
        
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
        end
    end
    
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = loglog(all_bin_centers{i}(which),pm_results.InputSNR(which,j),[':' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = loglog(all_bin_centers{i}(which),pm_results.InputSNR(which,j),[':' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
        % Determines which type of lines should be in the legend
        if j == 1 && i == numel(all_pm_results) 
            lines(end+1) = line;
        end
    end
    hold on;
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel('SNR (db)');
set(gca,'XScale','log');
legend(lines, clean_for_latex(legendentries2),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedSNRAxis')), ylim(TASBEConfig.get('OutputSettings.FixedSNRAxis')); end
title(['All ' clean_for_latex(stemName),' SNR vs. CFP']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-SNR'],directory);

% Delta SNR plots
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
lines = [];
for i=1:numel(all_pm_results)
    pm_results = all_pm_results{i};
    for j=1:step:n_var
        which = pm_results.Valid(:,j,1) & pm_results.Valid(:,j,2);
        if sum(which)==0, continue; end; % don't try to plot when lines are empty 
        if by_n_var
            line = semilogx(all_bin_centers{i}(which),pm_results.OutputSNR(which,j)-pm_results.InputSNR(which,j),['-' marker_types{i}],'Color',hsv2rgb([hues(j) 1 0.9])); hold on;
        else
            line = semilogx(all_bin_centers{i}(which),pm_results.OutputSNR(which,j)-pm_results.InputSNR(which,j),['-' marker_types{j}],'Color',hsv2rgb([hues(i) 1 0.9])); hold on;
        end
        % Determines which type of lines should be in the legend
        if j == 1 || i == 1 
            lines(end+1) = line;
        end
    end
    hold on;
end
xlabel(['CFP ' clean_for_latex(cfp_units)]); ylabel('\Delta SNR (db)');
set(gca,'XScale','log');
legend(lines, clean_for_latex(comlegendentries),'Location','Best');
if(TASBEConfig.isSet('OutputSettings.FixedBinningAxis')), xlim(TASBEConfig.get('OutputSettings.FixedBinningAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedDeltaSNRAxis')), ylim(TASBEConfig.get('OutputSettings.FixedDeltaSNRAxis')); end
title(['All ' clean_for_latex(stemName),' \Delta SNR vs. CFP']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-all-dSNR'],directory);
