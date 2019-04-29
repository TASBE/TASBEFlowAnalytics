% plot_autofluorescence_control generates a histogram showing the 
% how an AutoFluorescenceModel was generated from a non-fluorescent control
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plot_autofluorescence_control(AFM,data)

visiblePlots = TASBEConfig.get('autofluorescence.visiblePlots');
plotPath = TASBEConfig.get('autofluorescence.plotPath');
plotSize = TASBEConfig.get('autofluorescence.plotSize');

h = figure('PaperPosition',[1 1 plotSize]);
if(~visiblePlots), set(h,'visible','off'); end;
afmean = getMean(AFM);
afstd = getStd(AFM);
maxbin = max(500,afmean+2.1*afstd);
bins = BinSequence(-100, 10, maxbin, 'arithmetic');
bin_counts = zeros(size(get_bin_centers(bins)));
bin_edges = get_bin_edges(bins);
for k=1:numel(bin_counts)
    which = data>bin_edges(k) & data<=bin_edges(k+1);
    bin_counts(k) = sum(which);
end
plot(get_bin_centers(bins),bin_counts,'b-'); hold on;
plot([afmean afmean],[0 max(bin_counts)],'r-');
plot([afmean+2*afstd afmean+2*afstd],[0 max(bin_counts)],'r:');
plot([afmean-2*afstd afmean-2*afstd],[0 max(bin_counts)],'r:');
xlabel(sprintf('%s a.u.',clean_for_latex(getName(AFM.channel))));
ylabel('Count');
xlim([-100 maxbin]);
title(sprintf('Autofluorescence Model for %s',clean_for_latex(getPrintName(AFM.channel))));
outputfig(h, sprintf('autofluorescence-%s',getPrintName(AFM.channel)),plotPath);

