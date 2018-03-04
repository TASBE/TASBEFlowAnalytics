% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_inducer_characterization(results)

step = TASBEConfig.get('OutputSettings.PlotEveryN');
ticks = TASBEConfig.get('OutputSettings.PlotTickMarks');
stemName = TASBEConfig.get('OutputSettings.StemName');
deviceName = TASBEConfig.get('OutputSettings.DeviceName');
directory = TASBEConfig.get('plots.plotPath');

AP = getAnalysisParameters(results);
n_bins = get_n_bins(getBins(AP));
hues = (1:n_bins)./n_bins;

[input_mean input_std] = get_channel_results(results,'input');
in_units = getChannelUnits(AP,'input');

warning('TASBE:Plots','Assuming only a single inducer exists');
InducerName = getInducerName(getExperiment(results),1);
inducer_levels = getInducerLevelsToFiles(getExperiment(results),1);
which = inducer_levels==0;
if isempty(inducer_levels(inducer_levels>0)),
    inducer_levels(which) = 1;
else
    inducer_levels(which) = min(inducer_levels(inducer_levels>0))/10;
end
fa = getFractionActive(results);

%%%% Inducer plots
% Plain inducer plot:
h = figure('PaperPosition',[1 1 5 3.66]);
set(h,'visible','off');
for i=1:step:n_bins
    which = fa(i,:)>getMinFractionActive(AP);
    loglog(inducer_levels(which),input_mean(i,which),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    if ticks
        loglog(inducer_levels(which),input_mean(i,which),'+','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    end
    loglog(inducer_levels(which),input_mean(i,which).*input_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
    loglog(inducer_levels(which),input_mean(i,which)./input_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
end;
xlabel(['[',clean_for_latex(InducerName),']']); ylabel(['IFP ' clean_for_latex(in_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
if(TASBEConfig.isSet('OutputSettings.FixedInducerAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInducerAxis')); end;
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end;
title(['Raw ',clean_for_latex(deviceName),' transfer curve, colored by constitutive bin (non-equivalent colors)']);
outputfig(h,[clean_for_latex(stemName),'-',clean_for_latex(deviceName),'-mean'],directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasmid system is disabled, due to uncertainty about correctness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Normalized inducer plot:
% h = figure('PaperPosition',[1 1 5 3.66]);
% set(h,'visible','off');
% for i=1:step:n_bins
%     which = fa(i,:)>0.9;
%     pe=getPlasmidEstimates(results);
%     loglog(inducer_levels(which),input_mean(i,which)./pe(i,which),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
%     if ticks
%         loglog(inducer_levels(which),input_mean(i,which)./pe(i,which),'+','Color',hsv2rgb([hues(i) 1 0.9]));
%     end
%     loglog(inducer_levels(which),input_mean(i,which)./pe(i,which).*input_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
%     loglog(inducer_levels(which),input_mean(i,which)./pe(i,which)./input_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
% end;
% xlabel(['[',clean_for_latex(InducerName),']']); ylabel(['IFP ' clean_for_latex(in_units) '/plasmid']);
% set(gca,'XScale','log'); set(gca,'YScale','log');
% if(TASBEConfig.get('OutputSettings.FixedInducerAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInducerAxis')); end;
% if(TASBEConfig.get('OutputSettings.FixedNormalizedInputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedNormalizedInputAxis')); end;
% title(['Normalized ',TASBEConfig.get('OutputSettings.DeviceName'),' transfer curve, colored by plasmid bin (non-equivalent colors)']);
% outputfig(h,[TASBEConfig.get('OutputSettings.StemName'),'-',TASBEConfig.get('OutputSettings.DeviceName'),'-mean-norm'],TASBEConfig.get('OutputSettings.Directory'));
