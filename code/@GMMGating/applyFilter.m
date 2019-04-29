% APPLYFILTER overrides applyFilter method in Filter class. Filters raw FCS
% data for the GMMGating class. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function gated = applyFilter(GMMG,fcshdr,rawfcs)

n_channels = numel(GMMG.channel_names);

% gather channel data
channel_data = zeros(size(rawfcs,1),n_channels);
valid = ones(size(rawfcs,1),1);
for i=1:n_channels, 
    arith_channel_data = get_fcs_color(rawfcs,fcshdr,GMMG.channel_names{i});
    channel_data(:,i) = log10(arith_channel_data);
    valid = valid & arith_channel_data>0;
end;

% Compute clustering and distance matrix
clustered = cluster(GMMG.distribution,channel_data);
mh_dist_sq = mahal(GMMG.distribution,channel_data);

near = zeros(size(rawfcs,1),1);
for i=1:numel(GMMG.selected_components),
    batch = clustered==GMMG.selected_components(i);
    near = near | mh_dist_sq(:,GMMG.selected_components(i))<GMMG.deviations^2;
end

gated = rawfcs(valid & batch & near,:);
