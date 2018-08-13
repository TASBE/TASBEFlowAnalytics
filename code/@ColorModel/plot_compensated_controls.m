function plot_compensated_controls(CM)

% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

n = numel(CM.Channels);

if n==1, return; end; % nothing to plot if there's only one channel

for driven=1:n
    if(isUnprocessed(CM.Channels{driven})), continue; end; % skip unprocessed channels
    data = readfcs_compensated_au(CM,CM.ColorFiles{driven},0,0);

    for passive=1:n,
        if (passive == driven || isUnprocessed(CM.Channels{passive})), continue; end;
            
        h = figure('PaperPosition',[1 1 6 4]);
        set(h,'visible','off');
        pos = data(:,driven)>1;
        if sum(pos)==0, 
            TASBESession.warn('TASBE:CompensationModel','CannotCompensate','Cannot compensate %s with %s',getPrintName(CM.Channels{driven}),getPrintName(CM.Channels{passive}));
            continue; 
        end;
        pmin = percentile(data(pos,passive),0.1); pmax = percentile(data(pos,passive),99.9);
        smoothhist2D([log10(data(pos,driven)) data(pos,passive)],10,[200, 200],[],[],[0 pmin; 6 pmax]); hold on;
        %
        sorted = sortrows(data(pos,:),driven);
        for i=1:10, % passive means of deciles
            qsize = floor(size(sorted,1)/10);
            range = (1:qsize)+(qsize*(i-1));
            pmean = mean(sorted(range,passive));
            drange = [min(sorted(range,driven)) max(sorted(range,driven))];
            plot(log10(drange),[pmean pmean],'k*-');
        end
        
        cleanedDrivenPrintName = clean_for_latex(getPrintName(CM.Channels{driven}));
        cleanedPassivePrintName = clean_for_latex(getPrintName(CM.Channels{passive}));
        xlabel(sprintf('%s (%s a.u.)',cleanedDrivenPrintName,clean_for_latex(getName(CM.Channels{driven}))));
        ylabel(sprintf('%s (%s a.u.)',cleanedPassivePrintName,clean_for_latex(getName(CM.Channels{passive}))));
        title('Compensated Positive Control');
        path = TASBEConfig.get('compensation.plotPath');
        outputfig(h, sprintf('compensated-%s-vs-positive-%s',cleanedPassivePrintName,cleanedDrivenPrintName), path);
    end
end
