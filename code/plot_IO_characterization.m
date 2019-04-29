% PLOT_IO_CHARACTERIZATION plots plain input/output plots for transfer
% curve analysis.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_IO_characterization(results,in_channel,out_channel)
if nargin<2, in_channel = 'input'; end;
if nargin<3, out_channel = 'output'; end;

step = TASBEConfig.get('OutputSettings.PlotEveryN');
ticks = TASBEConfig.get('OutputSettings.PlotTickMarks');
stemName = TASBEConfig.get('OutputSettings.StemName');
directory = TASBEConfig.get('plots.plotPath');
figsize = TASBEConfig.get('OutputSettings.FigureSize');

AP = getAnalysisParameters(results);
n_bins = get_n_bins(getBins(AP));
hues = (1:n_bins)./n_bins;
fa = getFractionActive(results);

[input_mean] = get_channel_results(results,in_channel);
[output_mean output_std] = get_channel_results(results,out_channel);
in_units = getChannelUnits(AP,in_channel);
out_units = getChannelUnits(AP,out_channel);

%%% I/O plots:
% Plain I/O plot:
h = figure('PaperPosition',[1 1 figsize]);
set(h,'visible','off');
for i=1:step:n_bins
    which = fa(i,:)>getMinFractionActive(AP);
    loglog(input_mean(i,which),output_mean(i,which),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    if ticks
        loglog(input_mean(i,which),output_mean(i,which),'+','Color',hsv2rgb([hues(i) 1 0.9]));
    end
    % plot isolated points
    isolated = isolated_points(input_mean(i,:),1) | isolated_points(output_mean(i,:),1);
    loglog(input_mean(i,which & isolated),output_mean(i,which & isolated),'+','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
    % plot standard deviations
    loglog(input_mean(i,which),output_mean(i,which).*output_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
    loglog(input_mean(i,which),output_mean(i,which)./output_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
end
%if(TASBEConfig.get('OutputSettings.FixedAxis')), axis([1e2 1e10 1e2 1e10]); end;
xlabel(['IFP ' clean_for_latex(in_units)]); ylabel(['OFP ' clean_for_latex(out_units)]);
set(gca,'XScale','log'); set(gca,'YScale','log');
if(TASBEConfig.isSet('OutputSettings.FixedInputAxis')), xlim(TASBEConfig.get('OutputSettings.FixedInputAxis')); end
if(TASBEConfig.isSet('OutputSettings.FixedOutputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedOutputAxis')); end
title(['Raw ',clean_for_latex(stemName),' transfer curve, colored by constitutive bin']);
outputfig(h,[clean_for_latex(stemName),'-mean'],directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plasmid system is disabled, due to uncertainty about correctness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Normalized plot color by plasmid count:
% h = figure('PaperPosition',[1 1 5 3.66]);
% set(h,'visible','off');
% for i=1:step:n_bins
%     pe=getPlasmidEstimates(results);
%     which = fa(i,:)>0.9;
%     loglog(input_mean(i,which),output_mean(i,which)./pe(i,which),'-','Color',hsv2rgb([hues(i) 1 0.9])); hold on;
%     if ticks
%         loglog(input_mean(i,which),output_mean(i,which)./pe(i,which),'+','Color',hsv2rgb([hues(i) 1 0.9]));
%     end
%     loglog(input_mean(i,which),output_mean(i,which)./pe(i,which).*output_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
%     loglog(input_mean(i,which),output_mean(i,which)./pe(i,which)./output_std(i,which),':','Color',hsv2rgb([hues(i) 1 0.9]));
% end;
% %if(TASBEConfig.get('OutputSettings.FixedAxis')), axis([1e2 1e10 1e2 1e10]); end;
% set(gca,'XScale','log'); set(gca,'YScale','log');
% xlabel(['IFP ' clean_for_latex(in_units)]); ylabel(['OFP ' clean_for_latex(out_units) '/plasmid']);
% if(TASBEConfig.get('OutputSettings.FixedNormalizedInputAxis')), xlim(TASBEConfig.get('OutputSettings.FixedNormalizedInputAxis')); end;
% if(TASBEConfig.get('OutputSettings.FixedNormalizedOutputAxis')), ylim(TASBEConfig.get('OutputSettings.FixedNormalizedOutputAxis')); end;
% title(['Normalized ',TASBEConfig.get('OutputSettings.DeviceName'),' transfer curve, colored by plasmid count']);
% outputfig(h,[TASBEConfig.get('OutputSettings.StemName'),'-',TASBEConfig.get('OutputSettings.DeviceName'),'-mean-norm'],TASBEConfig.get('OutputSettings.Directory'));

end
