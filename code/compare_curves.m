% COMPARE_CURVES plots the difference between two characterization curves
% with the same inductions and segmented into the same subpopulation binning
%
% Curve 1 is shown as blue, Curve 2 as green
% Differences are shown as red lines between curves
% Points in only one data set are shown as black
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function compare_curves(in1,out1,in2,out2)
one_only = ~isnan(in1) & isnan(in2);
two_only = isnan(in1) & ~isnan(in2);
both = find(~isnan(in1) & ~isnan(in2));

mean_indiff = geomean(in1(both)./in2(both))
mean_outdiff = geomean(out1(both)./out2(both))

strongonly = find(~isnan(in1) & ~isnan(in2) & (in1>1e5 | 1e2>1e5));
mean_indiff_strong = geomean(in1(strongonly)./in2(strongonly))
mean_outdiff_strong = geomean(out1(strongonly)./out2(strongonly))

figsize = TASBEConfig.get('OutputSettings.FigureSize');
figure('PaperPosition',[1 1 figsize]);
loglog(in1(both),out1(both),'b.'); hold on;
loglog(in2(both),out2(both),'g.');
loglog(in1(one_only),out1(one_only),'k.');
loglog(in2(two_only),out2(two_only),'k.');
for i=1:numel(both)
    loglog([in1(both(i)) in2(both(i))],[out1(both(i)) out2(both(i))],'r-'); hold on;
end
xlabel('Input ERF'); ylabel('Output ERF');
TASBESession.warn('TASBE:Utilities','NotCheckingForPseudoERF','Not checking for pseudoERF');
