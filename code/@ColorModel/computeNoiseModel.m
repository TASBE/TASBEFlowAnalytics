% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function NM = computeNoiseModel(CM)

bins = BinSequence(3,0.1,9,'log_bins');

noisemins = zeros(numel(CM.Channels),1);
noisemeans = zeros(numel(CM.Channels),1);
noisestds = zeros(numel(CM.Channels),1);
detailcounts = cell(numel(CM.Channels),1);
detailmeans = cell(numel(CM.Channels),1);
detailstds = cell(numel(CM.Channels),1);

for i=1:numel(CM.Channels)
    clear counts means stds minset;
    n_files = 0;
    
    for j=1:numel(CM.ColorPairFiles)
        cp = CM.ColorPairFiles{j};
        % Find the focus channel in the pair-files
        found = false; target_id = 3;
        for k=1:2 % find target channel
            if(CM.Channels{i}==cp{k})
                target = indexof(CM.Channels,cp{k}); found=true; target_id = k;
                if(target==-1), error('Missing channel %s',getPrintName(cp{k})); end
            end
        end
        if ~found, continue; end;
        ctrl_idx = 3;
        if(cp{3}==cp{target_id}) % two-color case
            ctrl_idx = 3-target_id;
        end
        cCtrl = indexof(CM.Channels,cp{ctrl_idx});
        if(cCtrl==-1), error('Missing channel %s',getPrintName(cp{ctrl_idx})); end
        
        % read the file and get another layer of statistics
        data = readfcs_compensated_ERF(CM,cp{4},false,true);
        [localcounts localmeans localstds] = subpopulation_statistics(bins,data,cCtrl,'geometric');
        n_files=n_files+1;
        counts(:,n_files) = localcounts;
        means(:,n_files) = localmeans(:,target);
        stds(:,n_files) = localstds(:,target);
        detailcounts{i}{n_files} = counts(:,n_files);
        detailmeans{i}{n_files} = means(:,n_files);
        detailstds{i}{n_files} = stds(:,n_files);

        which = localmeans(:,target)>(getStdERF(CM.autofluorescence_model{i})*2) & localcounts>100;
        if(sum(which)>0)
            minset(n_files) = min(localstds(which,target));
        else
            warning('Model:Noise','Color mapping file "%s" has no significant fluorescent subpopulations',cp{4});
            minset(n_files) = 0;
        end
    end
    
    if n_files == 0,
        warning('Model:Noise','Cannot compute noise model for %s: no color mapping files available',getName(CM.Channels{i}));
        continue;
    end
    
    which = (means>getStdERF(CM.autofluorescence_model{i})*2) & counts>100;
    noisemins(i) = mean(minset);
    noisemeans(i) = mean(stds(which));
    noisestds(i) = std(stds(which));
    
    if CM.noise_plot
        h = figure('PaperPosition',[1 1 6 4]);
        set(h,'visible','off');
        barmin = 1e9; barmax = 1e3;
        for j=1:n_files
            which = means(:,j)>getStdERF(CM.autofluorescence_model{i})*2 & counts(:,j)>100;
            semilogx(means(which,j),stds(which,j),'b-'); hold on;
            barmin = min(min(means(which,j)),barmin); barmax = max(max(means(which,j)),barmax);
        end
        plot([1e3 1e9],[noisemins(i) noisemins(i)],'g-');
        plot([1e3 1e9],[noisemeans(i) noisemeans(i)],'r-');
        plot([1e3 1e9],[noisemeans(i)+noisestds(i) noisemeans(i)+noisestds(i)],'r:');
        plot([1e3 1e9],[noisemeans(i)-noisestds(i) noisemeans(i)-noisestds(i)],'r:');
        xlabel(sprintf('%s %s',clean_for_latex(getName(CM.Channels{i})),clean_for_latex(getUnits(CM.Channels{i}))));
        ylabel('k-fold noise');
        title(sprintf('Noise Model for %s',clean_for_latex(getName(CM.Channels{i}))));
        path = TASBEConfig.get('plots.plotPath');
        outputfig(h,sprintf('noise-model-%s',clean_for_latex(getPrintName(CM.Channels{i}))),path);
    end
end

NM = NoiseModel(noisemins,noisemeans,noisestds,detailcounts,detailmeans,detailstds);
