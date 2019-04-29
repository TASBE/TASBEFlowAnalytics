% function GMMG = GMMGating(file)
%   Constructor of GMMGating class, which is a subclass of Filter
%   file: a datafile or string for the data to be used for gating
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function GMMG = GMMGating(file)

GMMG.selected_components = [];
GMMG.channel_names = {};
GMMG.distribution = {};
GMMG.deviations = [];
GMMG.fraction_kept = 0.0; % will report the fraction of data that weren't filtered out'

% gate function just runs autogate_filter on model
GMMG = class(GMMG,'GMMGating',Filter());

if nargin==0, return; end;

% Model is a gmdistribution

file = ensureDataFile(file);

[~, fcshdr, rawfcs] = fca_read(file);

channel_names = TASBEConfig.get('gating.channelNames');
n_channels = numel(channel_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain and prefilter data

gate_fraction = TASBEConfig.get('gating.fraction');
k_components = TASBEConfig.get('gating.kComponents');
selected_components = TASBEConfig.getexact('gating.selectedComponents',[]);
gate_deviations = TASBEConfig.get('gating.deviations');
gate_tightening = TASBEConfig.get('gating.tightening');

% gather channel data
unfiltered_channel_data = cell(n_channels,1);
unfiltered_channel_data_arith = unfiltered_channel_data;
for i=1:n_channels 
    unfiltered_channel_data_arith{i} = get_fcs_color(rawfcs,fcshdr,channel_names{i});
    unfiltered_channel_data{i} = log10(unfiltered_channel_data_arith{i});
end

% filter channel data away from saturation points
which = ones(numel(unfiltered_channel_data{1}),1);
for i=1:n_channels
    valid = ~isinf(unfiltered_channel_data{i}) & ~isnan(unfiltered_channel_data{i}) & (unfiltered_channel_data_arith{i}>0);
    bound = [min(unfiltered_channel_data{i}(valid)) max(unfiltered_channel_data{i}(valid))];
    span = bound(2)-bound(1);
    range = [mean(bound)-span*gate_fraction/2 mean(bound)+span*gate_fraction/2];
    which = which & valid & unfiltered_channel_data{i}>range(1) & unfiltered_channel_data{i}<range(2);
end
channel_data = zeros(sum(which),n_channels);
for i=1:n_channels 
    channel_data(:,i) = unfiltered_channel_data{i}(which);
end
frac_kept = sum(which)/numel(which);
fprintf('Gating autodetect using %.2f%% valid and non-saturated data\n',100*frac_kept);
GMMG.fraction_kept = frac_kept;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find and adjust gaussian fit
% Control random seed if fixedSeed TASBEConfig is true
if TASBEConfig.get('gating.fixedSeed')
    if is_octave
        rand ('seed', 10);
        dist = fitgmdist(channel_data,k_components,'Regularize',1e-5);
    else
        rng(10); % For reproducibility
        dist = fitgmdist(channel_data,k_components,'Regularize',1e-5);
    end
else
    % If fixedSeed is false, call fitgmdist as normal
    dist = fitgmdist(channel_data,k_components,'Regularize',1e-5);
end

% assignin('base','GMMdist',dist);
dss = struct(dist); %% Terrible kludge: should actually make accessors
% sort component identities by eigenvalue size
maxeigs = zeros(k_components,1);
for i=1:k_components
    maxeigs(i) = max(eig(dss.Sigma(:,:,i)));
end
sorted_eigs = sortrows([maxeigs'; 1:k_components]');
eigsort = sorted_eigs(:,2);

if isempty(selected_components)
    selected_components = 1;
end
GMMG.selected_components = eigsort(selected_components);

% reweight components to make select components tighter
reweight = dss.PComponents;
lossweight = gate_tightening*sum(reweight(GMMG.selected_components));
for i=1:k_components,
    if(isempty(find(i==GMMG.selected_components, 1)))
        reweight(i) = reweight(i)-gate_tightening*reweight(i);
    else
        reweight(i) = reweight(i)*(1+lossweight);
    end
end

% Assembly GMMG package:
GMMG.channel_names = channel_names;
GMMG.distribution = gmdistribution(dss.mu,dss.Sigma,reweight);
GMMG.deviations = gate_deviations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the plots
%%%%
makePlots = TASBEConfig.get('gating.plot');
visiblePlots = TASBEConfig.get('gating.visiblePlots');
plotPath = TASBEConfig.get('gating.plotPath');
plotSize = TASBEConfig.get('gating.plotSize');
showNonselected = TASBEConfig.get('gating.showNonselected');
largeOutliers = TASBEConfig.get('gating.largeOutliers');
range = TASBEConfig.getexact('gating.range',[]);
density = TASBEConfig.get('gating.density');

if makePlots
    if density >= 1, type = 'image'; else type = 'contour'; end

    % Plot vs. gate:
    gated = applyFilter(GMMG,fcshdr,rawfcs);
    % compute component gates too
    if showNonselected,
        fprintf('Computing individual components');
        c_gated = cell(k_components,1);
        for i=1:k_components
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
        h = figure('PaperPosition',[1 1 plotSize]);
        if(~visiblePlots), set(h,'visible','off'); end;
        %smoothhist2D([channel_data{1} channel_data{2}],10,[200, 200],[],type,range,largeOutliers);
        smoothhist2D([channel_data(:,i) channel_data(:,i+1)],5,[500, 500],[],type,range,largeOutliers);
        xlabel([clean_for_latex(channel_names{i}) ' a.u.']); 
        ylabel([clean_for_latex(channel_names{i+1}) ' a.u.']);
        title('2D Gaussian Gate Fit');
        hold on;

        % Show components:
        if showNonselected,
            component_color = [0.5 0.0 0.0];
            component_text = [0.8 0.0 0.0];
            for j=1:k_components
                gated_sub = [log10(get_fcs_color(c_gated{j},fcshdr,channel_names{i})) log10(get_fcs_color(c_gated{j},fcshdr,channel_names{i+1}))];
                if size(gated_sub,1) > 3,
                    which = convhull(gated_sub(:,1),gated_sub(:,2));
                    plot(gated_sub(which,1),gated_sub(which,2),'-','LineWidth',2,'Color',component_color);
                    text_x = mean(gated_sub(which,1)); text_y = mean(gated_sub(which,2));
                    text(text_x,text_y,sprintf('Component %i',j),'Color',component_text);
                end
            end
        end

        % Show gated:
        gated_sub = [log10(get_fcs_color(gated,fcshdr,channel_names{i})) log10(get_fcs_color(gated,fcshdr,channel_names{i+1}))];
        which = convhull(gated_sub(:,1),gated_sub(:,2));
        plot(gated_sub(which,1),gated_sub(which,2),'r-','LineWidth',2);

        outputfig(h,clean_for_latex(sprintf('AutomaticGate-%s-vs-%s',channel_names{i},channel_names{i+1})), plotPath);
    end
end
