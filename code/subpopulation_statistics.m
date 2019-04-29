% SUBPOPULATION_STATISTICS returns the stats (counts, means, stds, and excluded) of a certain subset of
% inputted data determined by the selector input. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [counts, means, stds, excluded] = subpopulation_statistics(BSeq,data,selector,mode,drop_thresholds)
if nargin<5, drop_thresholds = []; end;

bedges = get_bin_edges(BSeq);

n = numel(bedges)-1;
ncol = size(data,2);

nbins = get_n_bins(BSeq);

% create zero sets
counts = zeros(nbins,1); 
means = zeros(nbins,ncol); stds = means;

switch(mode)
    case 'geometric'
        for i=1:n
            selection = data(:,selector)>bedges(i) & data(:,selector)<=bedges(i+1);
            if ~isempty(drop_thresholds), selection = selection & data(:,selector)>drop_thresholds(selector); end;
            which = find(selection);
            counts(i) = numel(which);

            if counts(i)==0,
                means(i,:) = NaN;
                stds(i,:) = NaN;
            else
                for j=1:ncol
                    trimmed = data(which,j);

                    means(i,j) = geomean(trimmed);
                    stds(i,j) = geostd(trimmed);
                end
            end
        end
      
    case 'arithmetic'
        for i=1:n
            which = find(data(:,selector)>bedges(i) & data(:,selector)<=bedges(i+1));
            counts(i) = numel(which);
            
            if counts(i)==0,
                means(i,:) = NaN;
                stds(i,:) = NaN;
            else
                for j=1:ncol
                    trimmed = data(which,j);

                    means(i,j) = mean(trimmed);
                    stds(i,j) = std(trimmed);
                end
            end
        end
        
    otherwise
        TASBESession.error('TASBE:Analysis','UnknownStatisticalMode','Unknown statistical mode %s',mode);
end

excluded = size(data,1) - sum(counts);

end
