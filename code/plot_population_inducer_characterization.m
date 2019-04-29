% PLOT_POPULATION_INDUCER_CHARACTERIZATION creates plain inducer plots for transfer curve analysis. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_population_inducer_characterization(results)

ticks = TASBEConfig.get('OutputSettings.PlotTickMarks');
stemName = TASBEConfig.get('OutputSettings.StemName');
directory = TASBEConfig.get('plots.plotPath');
deviceName = TASBEConfig.get('OutputSettings.DeviceName');
figsize = TASBEConfig.get('OutputSettings.FigureSize');

AP = getAnalysisParameters(results);
n_components = getNumGaussianComponents(AP);
hues = (1:n_components)./n_components;

[input_mean input_std] = get_channel_population_results(results,'input');
in_units = getChannelUnits(AP,'input');

TASBESession.warn('TASBE:Plots','SingleInducerOnly','Assuming only a single inducer exists');
InducerName = getInducerName(getExperiment(results),1);
inducer_levels = getInducerLevelsToFiles(getExperiment(results),1);
which = inducer_levels==0;
% Find the smallest non-zero value, min value, and max_value and thresholds from inducer_levels
min_non_zero = min(inducer_levels(inducer_levels>0));
min_value = min(inducer_levels);
max_value = max(inducer_levels);
higher_threshold = log10(min_non_zero)-1; % stands for start in ZeroOnLog function
lower_threshold = higher_threshold -1; % stands for zero in ZeroOnLog function
% set the zero values to the lower_threshold
inducer_levels(which) = 10^lower_threshold;

%%%% Inducer plots
% Plain inducer plot:
h = figure('PaperPosition',[1 1 figsize]);
set(h,'visible','off');
for i=1:n_components
    loglog(inducer_levels(:),10.^input_mean(i,:),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    if ticks
        loglog(inducer_levels(:),10.^input_mean(i,:),'+','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    end
    loglog(inducer_levels(:),10.^(input_mean(i,:)+input_std(i,:)),':','Color',hsv2rgb([hues(i) 1 0.9]));
    loglog(inducer_levels(:),10.^(input_mean(i,:)-input_std(i,:)),':','Color',hsv2rgb([hues(i) 1 0.9]));
end;
xlabel(['[',clean_for_latex(InducerName),']']); ylabel(['IFP ' clean_for_latex(in_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
if(TASBEConfig.isSet('OutputSettings.FixedInducerAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInducerAxis')); end;
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end;
title(['Population ',clean_for_latex(deviceName),' transfer curve, colored by Gaussian component']);
% Edit ticks on plot to include 0 
if min_value <= 0
    xlim([10^lower_threshold max_value]); % set limit of x-axis to match with zero and start positions
    ZeroOnLog(10^lower_threshold,0.5*10^higher_threshold); % call ZeroOnLog with thresholds (0.5 is for scaling of '\\')
end
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-pop-mean'],directory);
