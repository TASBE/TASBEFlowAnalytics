% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function autofluorescence_model=computeAutoFluorescence(CM)

path = TASBEConfig.get('plots.plotPath');

% Import data from FCS files
[rawfcs fcshdr] = read_filtered_au(CM,CM.BlankFile);

autofluorescence_model = cell(numel(CM.Channels),1);
for i=1:numel(CM.Channels)
    found = false;
    name=getName(CM.Channels{i});
    for j=1:numel(fcshdr.par)
        if(strcmp(name,fcshdr.par(j).name) || strcmp(name,fcshdr.par(j).rawname))
            autofluorescence_model{i} = AutoFluorescenceModel(rawfcs(:,j));
            if(CM.Channels{i} == CM.ERF_channel)
                autofluorescence_model{i}=ERFize(autofluorescence_model{i},1,getK_ERF(CM.unit_translation));
            end
            
            % Optional plot
            if CM.autofluorescence_plot
                h = figure('PaperPosition',[1 1 6 4]);
                set(h,'visible','off');
                afmean = getMean(autofluorescence_model{i});
                afstd = getStd(autofluorescence_model{i});
                maxbin = max(500,afmean+2.1*afstd);
                bins = BinSequence(-100, 10, maxbin, 'arithmetic');
                bin_counts = zeros(size(get_bin_centers(bins)));
                bin_edges = get_bin_edges(bins);
                for k=1:numel(bin_counts)
                    which = rawfcs(:,j)>bin_edges(k) & rawfcs(:,j)<=bin_edges(k+1);
                    bin_counts(k) = sum(which);
                end
                plot(get_bin_centers(bins),bin_counts,'b-'); hold on;
                plot([afmean afmean],[0 max(bin_counts)],'r-');
                plot([afmean+2*afstd afmean+2*afstd],[0 max(bin_counts)],'r:');
                plot([afmean-2*afstd afmean-2*afstd],[0 max(bin_counts)],'r:');
                xlabel(sprintf('%s a.u.',clean_for_latex(name)));
                ylabel('Count');
                xlim([-100 maxbin]);
                title('Autofluorescence Model');
                outputfig(h, sprintf('autofluorescence-%s',clean_for_latex(getPrintName(CM.Channels{i}))),path);
            end
            
            found = true; break;
        end
    end
    if(found==false)
        error('Unable to find required channel %s',name)
    end
end
