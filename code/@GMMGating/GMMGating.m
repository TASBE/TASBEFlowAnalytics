% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function GMMG = GMMGating(file, AGP, output_path)

GMMG.selected_components = [];
GMMG.channel_names = {};
GMMG.distribution = {};
GMMG.deviations = [];
GMMG.fraction_kept = 0.0; % will report the fraction of data that weren't filtered out'

% gate function just runs autogate_filter on model
GMMG = class(GMMG,'GMMGating',Filter());

if nargin==0, return; end;

% AGP = AutogateParameters
% Model is a gmdistribution

[~, fcshdr, rawfcs] = fca_readfcs(file);

if(nargin<2), AGP = AutogateParameters(); end;
if (nargin < 3)
    output_path = './';
end


n_channels = numel(AGP.channel_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain and prefilter data

% gather channel data
unfiltered_channel_data = cell(n_channels,1);
unfiltered_channel_data_arith = unfiltered_channel_data;
for i=1:n_channels, 
    unfiltered_channel_data_arith{i} = get_fcs_color(rawfcs,fcshdr,AGP.channel_names{i});
    unfiltered_channel_data{i} = log10(unfiltered_channel_data_arith{i});
end;

% filter channel data away from saturation points
which = ones(numel(unfiltered_channel_data{1}),1);
for i=1:n_channels,
    valid = ~isinf(unfiltered_channel_data{i}) & ~isnan(unfiltered_channel_data{i}) & (unfiltered_channel_data_arith{i}>0);
    bound = [min(unfiltered_channel_data{i}(valid)) max(unfiltered_channel_data{i}(valid))];
    span = bound(2)-bound(1);
    range = [mean(bound)-span*AGP.fraction/2 mean(bound)+span*AGP.fraction/2];
    which = which & valid & unfiltered_channel_data{i}>range(1) & unfiltered_channel_data{i}<range(2);
end
channel_data = zeros(sum(which),n_channels);
for i=1:n_channels, 
    channel_data(:,i) = unfiltered_channel_data{i}(which);
end;
frac_kept = sum(which)/numel(which);
fprintf('Gating autodetect using %.2f%% valid and non-saturated data\n',100*frac_kept);
GMMG.fraction_kept = frac_kept;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and adjust gaussian fit

dist = fitgmdist(channel_data,AGP.k_components,'Regularize',1e-5);
dss = struct(dist); %% Terrible kludge: should actually make accessors
% sort component identities by eigenvalue size
maxeigs = zeros(AGP.k_components,1);
for i=1:AGP.k_components,
    maxeigs(i) = max(eig(dss.Sigma(:,:,i)));
end
sorted_eigs = sortrows([maxeigs'; 1:AGP.k_components]');
eigsort = sorted_eigs(:,2);

if isempty(AGP.selected_components),
    AGP.selected_components = 1;
end
GMMG.selected_components = eigsort(AGP.selected_components);

% reweight components to make select components tighter
reweight = dss.PComponents;
lossweight = AGP.tightening*sum(reweight(GMMG.selected_components));
for i=1:AGP.k_components,
    if(isempty(find(i==GMMG.selected_components, 1)))
        reweight(i) = reweight(i)-AGP.tightening*reweight(i);
    else
        reweight(i) = reweight(i)*(1+lossweight);
    end
end

% Assembly GMMG package:
GMMG.channel_names = AGP.channel_names;
GMMG.distribution = gmdistribution(dss.mu,dss.Sigma,reweight);
GMMG.deviations = AGP.deviations;

% AGP = AutogateParameters
% Model is a gmdistribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the plots
%%%%
if AGP.density >= 1, type = 'image'; else type = 'contour'; end

% Plot vs. gate:
gated = applyFilter(GMMG,fcshdr,rawfcs);
% compute component gates too
if AGP.show_nonselected,
    fprintf('Computing individual components');
    c_gated = cell(AGP.k_components,1);
    for i=1:AGP.k_components
        tmp_model = GMMG;
        tmp_model.selected_components = eigsort(i);
        c_gated{i} = applyFilter(tmp_model,fcshdr,rawfcs);
        fprintf('.');
    end
    fprintf('\n');
end

for i=1:2:n_channels,
    % handle odd number of channels by decrementing last:
    if i==n_channels, i=i-1; end;
    
    % Show background:
    h = figure('PaperPosition',[1 1 5 5]);
    if ~AGP.visible, set(h,'visible','off'); end;
    %smoothhist2D([channel_data{1} channel_data{2}],10,[200, 200],[],type,range,largeoutliers);
    smoothhist2D([channel_data(:,i) channel_data(:,i+1)],5,[500, 500],[],type,AGP.range,AGP.largeoutliers);
    xlabel(AGP.channel_names{i}); ylabel(AGP.channel_names{i+1});
    title('2D Gaussian Gate Fit');
    hold on;

    % Show components:
    if AGP.show_nonselected,
        component_color = [0.5 0.0 0.0];
        component_text = [0.8 0.0 0.0];
        for j=1:AGP.k_components
            gated_sub = [log10(get_fcs_color(c_gated{j},fcshdr,AGP.channel_names{i})) log10(get_fcs_color(c_gated{j},fcshdr,AGP.channel_names{i+1}))];
            if size(gated_sub,1) > 3,
                which = convhull(gated_sub(:,1),gated_sub(:,2));
                plot(gated_sub(which,1),gated_sub(which,2),'-','LineWidth',2,'Color',component_color);
                text_x = mean(gated_sub(which,1)); text_y = mean(gated_sub(which,2));
                text(text_x,text_y,sprintf('Component %i',j),'Color',component_text);
            end
        end
    end
    
    % Show gated:
    gated_sub = [log10(get_fcs_color(gated,fcshdr,AGP.channel_names{i})) log10(get_fcs_color(gated,fcshdr,AGP.channel_names{i+1}))];
    which = convhull(gated_sub(:,1),gated_sub(:,2));
    plot(gated_sub(which,1),gated_sub(which,2),'r-','LineWidth',2);

    outputfig(h,sprintf('AutomaticGate-%s-vs-%s',AGP.channel_names{i},AGP.channel_names{i+1}), output_path);
end
