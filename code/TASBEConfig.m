% Copyright (C) 2011 - 2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE Flow Analytics distribution's top directory.
%
% This file is part of the TASBE Flow Analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the BBN Flow Cytometry
% package distribution's top directory.

classdef TASBEConfig
    methods(Static,Hidden)
        function state = init()
            s = struct();
            doc = struct();
            doc.about = 'Configuration variables for TASBE Flow Analytics';
            defaults = containers.Map();
            
            % Testing variables
            s.testing = struct(); doc.testing = struct();
            doc.testing.about = 'Settings to disable features for testing purposes';
            doc.testing.fakeFigureSaves = 'For testing purposes, do not actually save figures, only pretend to.';
            s.testing.fakeFigureSaves = 0;
            
            % Matlab GMdistribution
            
            % Generic flow data analysis
            s.flow = struct(); doc.flow = struct();
            doc.flow.about = 'General settings for flow cytometry data analysis';
            doc.flow.maxEvents = 'Drop events above this count to avoid memory issues';
            s.flow.maxEvents = 1e6;
            doc.flow.rangeMin = 'bin minimum (log10 scale)';
            s.flow.rangeMin = 0;                           
            doc.flow.rangeMax = 'bin maximum (log10 scale)';
            s.flow.rangeMax = 7;
            doc.flow.onThreshold = 'Threshold to differentiate between on and off events. Points with values equal to threshold are placed in the on group.';
            s.flow.onThreshold = 0;
            doc.flow.outputPointCloud = 'if true, output point-cloud for each calibrated read';
            s.flow.outputPointCloud = false;
            doc.flow.pointCloudPath = 'location for point-cloud outputs';
            s.flow.pointCloudPath = 'CSV/';
            doc.flow.pointCloudFileType = 'what type of pathname to write to point-cloud JSON header file. 0 stands for absolute, 1 stands for relative.';
            s.flow.pointCloudFileType = 0;
            doc.flow.dataCSVPath = 'location for data summary CSVs';
            s.flow.dataCSVPath = 'CSV/';
            doc.flow.outputHistogramFile = 'if true, output histogram file for batch analysis';
            s.flow.outputHistogramFile = true;
            doc.flow.outputStatisticsFile = 'if true, output statistics file for batch analysis';
            s.flow.outputStatisticsFile = true;
            doc.flow.channel_template_file = 'File to check laser/filter settings against (defaults to bead file)';
            s.flow.channel_template_file = [];
            % sample size/gating warning thresholds
            doc.flow.smallFileWarning = 'Threshold for warning that event count is unusually low';
            s.flow.smallFileWarning = 10000;
            doc.flow.gateDiscardsWarning = 'Threshold for warning gates are discarding too high a fraction of events';
            s.flow.gateDiscardsWarning = 0.1;
            doc.flow.preGateDiscardsWarning = 'Threshold for warning pre-gates are discarding too high a fraction of events';
            s.flow.preGateDiscardsWarning = [];
            defaults('flow.preGateDiscardsWarning') = 'flow.gateDiscardsWarning';
            doc.flow.postGateDiscardsWarning = 'Threshold for warning post-gates are discarding too high a fraction of events';
            s.flow.postGateDiscardsWarning = [];
            defaults('flow.postGateDiscardsWarning') = 'flow.gateDiscardsWarning';
            doc.flow.eventRatioWarning = 'Generic threshold for warning about different event counts.';
            s.flow.eventRatioWarning = 10;
            doc.flow.replicateEventRatioWarning = 'Threshold to warn on variation in replicate event counts.';
            s.flow.replicateEventRatioWarning = [];
            defaults('flow.replicateEventRatioWarning') = 'flow.eventRatioWarning';
            doc.flow.conditionEventRatioWarning = 'Threshold to warn on variation in condition event counts.';
            s.flow.conditionEventRatioWarning = [];
            defaults('flow.conditionEventRatioWarning') = 'flow.eventRatioWarning';

            % generic plots
            s.plots = struct(); doc.plots = struct();
            doc.plots.about = 'General settings for plotting figures';
            doc.plots.plotPath = 'Default location for plots';
            s.plots.plotPath = 'plots/';
            doc.plots.visiblePlots = 'If true, plots are visible; otherwise, they are hidden for later saving';
            s.plots.visiblePlots = false;
            doc.plots.graphPlotSize = 'Default size (in inches) [X Y] for data graph figures';
            s.plots.graphPlotSize = [6 4];
            doc.plots.heatmapPlotSize = 'Default size (in inches) [X Y] for scatter/heatmap figures';
            s.plots.heatmapPlotSize = [5 5];
            doc.plots.largeOutliers = 'If true, outliers in heatmap figures are large, for output in small figures';
            s.plots.largeOutliers = false;
            doc.plots.heatmapPlotType = 'Set to ''image'', ''contour'', or ''surf'' to determine type of heatmap images';
            s.plots.heatmapPlotType = 'image';
            
            % calibration plots, i.e., those supporting the transformation of raw data into processed data, like autofluorescence, compensation, units
            s.calibration = struct(); doc.calibration = struct();
            doc.calibration.about = 'General settings for calibration figures';
            doc.calibration.overrideUnits = 'When set, force a.u. to ERF scaling value to this value';
            s.calibration.overrideUnits = [];
            doc.calibration.overrideAutofluorescence = 'When set, force autofluorescence to use this as mean AF';
            s.calibration.overrideAutofluorescence = [];
            doc.calibration.overrideCompensation = 'When set, use this matrix instead of computing a linear compensation model';
            s.calibration.overrideCompensation = [];
            doc.calibration.overrideTranslation = 'When set, use this matrix instead of computing a color translation model';
            s.calibration.overrideTranslation = [];
            doc.calibration.plot = 'When true, make diagnostic plots while computing color models';
            s.calibration.plot = true;
            doc.calibration.visiblePlots = 'If true, calibration plots are visible; otherwise, they are hidden for later saving';
            s.calibration.visiblePlots = [];
            defaults('calibration.visiblePlots') = 'plots.visiblePlots';
            doc.calibration.plotPath = 'Default location for calibration plots';
            s.calibration.plotPath = []; 
            defaults('calibration.plotPath') = 'plots.plotPath';
            doc.calibration.graphPlotSize = 'Default size (in inches) [X Y] for calibration data graph figures';
            s.calibration.graphPlotSize = [];
            defaults('calibration.graphPlotSize') = 'plots.graphPlotSize';
            doc.calibration.heatmapPlotSize = 'Default size (in inches) [X Y] for calibration scatter/heatmap figures';
            s.calibration.heatmapPlotSize = [];
            defaults('calibration.heatmapPlotSize') = 'plots.heatmapPlotSize';
            
            % Gating
            s.gating = struct(); doc.gating = struct();
            doc.gating.about = 'Settings for GMM Gating';
            doc.gating.fixedSeed = 'When true, controls the random seed for GMM Gating';
            s.gating.fixedSeed = true;
            % gating control parameters
            doc.gating.deviations = 'Number of standard deviations within which gate selects data';
            s.gating.deviations = 2.0;
            doc.gating.tightening = 'Amount that selected gate components are further tightened (range: [0,1])';
            s.gating.tightening = 0.0;
            doc.gating.kComponents = 'Number of gaussian components in a gate';
            s.gating.kComponents = 2;
            doc.gating.channelNames = 'Channels to use for gating'; % used to also include 'FSC-H','FSC-W','SSC-H','SSC-W'
            s.gating.channelNames = {'FSC-A','SSC-A'}; 
%             s.gating.saturationWarning = 0.3;       % Warn about saturation distorting gates if fraction non-extrema less than this level
            doc.gating.fraction = 'Fraction of range considered saturated and thus excluded from computation';
            s.gating.fraction = 0.95;
            doc.gating.selectedComponents = 'Set to force selection of particular components';
            s.gating.selectedComponents = [];
            % gating custom plotting parameters
            doc.gating.showNonselected = 'When true, show all gate components; when false, selected components only.';
            s.gating.showNonselected = true;
            doc.gating.largeOutliers = 'When true, plot gate outliers with large dots';
            s.gating.largeOutliers = true;
            doc.gating.range = 'Force gate heatmap plotting range of [x_min y_min; x_max y_max]';
            s.gating.range = [];
            doc.gating.density = 'Gate heatmap style: set to 1 for image, 0 for contour';
            s.gating.density = 1;
            % gating standard plotting parameters
            doc.gating.plot = 'Determines whether gating plots should be created';
            s.gating.plot = [];
            defaults('gating.plot') = 'calibration.plot';
            doc.gating.visiblePlots = 'If true, gating plots are visible; otherwise, they are hidden for later saving';
            s.gating.visiblePlots = [];
            defaults('gating.visiblePlots') = 'calibration.visiblePlots';
            doc.gating.plotPath = 'Location for gating plots';
            s.gating.plotPath = [];                 % where should gating plot go?
            defaults('gating.plotPath') = 'calibration.plotPath';
            doc.gating.plotSize = 'Size (in inches) [X Y] for gating figures';
            s.gating.plotSize = [];                 % What size (in inches) should gating plot be?
            defaults('gating.plotSize') = 'calibration.heatmapPlotSize';
            
            
            % Autofluorescence removal
            s.autofluorescence = struct(); doc.autofluorescence = struct();
            doc.autofluorescence.about = 'Settings for autofluorescence models';
            doc.autofluorescence.dropFraction = 'Fraction of extrema (high and low) to drop before computing autofluorescence';
            s.autofluorescence.dropFraction = 0.025;
            doc.autofluorescence.plot = 'Determines whether autofluorescence plots should be created';
            s.autofluorescence.plot = [];
            defaults('autofluorescence.plot') = 'calibration.plot';
            doc.autofluorescence.visiblePlots = 'If true, autofluorescence plots are visible; otherwise, they are hidden for later saving';
            s.autofluorescence.visiblePlots = [];
            defaults('autofluorescence.visiblePlots') = 'calibration.visiblePlots';
            doc.autofluorescence.plotPath = 'Location for autofluorescence plots';
            s.autofluorescence.plotPath = [];
            defaults('autofluorescence.plotPath') = 'calibration.plotPath';
            doc.autofluorescence.plotSize = 'Size (in inches) [X Y] for autofluorescence figures';
            s.autofluorescence.plotSize = [];
            defaults('autofluorescence.plotSize') = 'calibration.graphPlotSize';
            
            % Spectral bleed compensation
            s.compensation = struct(); doc.compensation = struct();
            doc.compensation.about = 'Settings for spectral bleed compensation';
            doc.compensation.minimumDrivenLevel = 'Uniformly ignores all less than this level of a.u.';
            s.compensation.minimumDrivenLevel = 1e2;    % uniformly ignore all less than this level of a.u. 
            doc.compensation.maximumDrivenLevel = 'Uniformly ignores all greater than this level of a.u.';
            s.compensation.maximumDrivenLevel = Inf;     % uniformly ignore all greater than this level of a.u. 
            doc.compensation.minimumBinCount = 'Ignores bins with less than this many elements';
            s.compensation.minimumBinCount = 10;        % ignore bins with less than this many elements
            doc.compensation.highBleedWarning = 'Warns about high bleed at this level';
            s.compensation.highBleedWarning = 0.1;      % Warn about high bleed at this level
            doc.compensation.plot = 'Determines whether compensation plots should be created';
            s.compensation.plot = [];                   % Should compensation plots be created?
            defaults('compensation.plot') = 'calibration.plot';
            doc.compensation.visiblePlots = 'If true, compensation plots are visible; otherwise, they are hidden for later saving';
            s.compensation.visiblePlots = [];           % should compensation plot be visible, or just created?
            defaults('compensation.visiblePlots') = 'calibration.visiblePlots';
            doc.compensation.plotPath = 'Location for compensation plots';
            s.compensation.plotPath = [];               % where should compensation plot go?
            defaults('compensation.plotPath') = 'calibration.plotPath';
            doc.compensation.plotSize = 'Size (in inches) [X Y] for compensation figures';
            s.compensation.plotSize = [6 4];
            defaults('compensation.plotSize') = 'calibration.heatmapPlotSize';
            
            % Color translation
            s.colortranslation = struct();
            doc.colortranslation.about = 'Settings for color translation models';
            doc.colortranslation.rangeMin = 'Minimum for color translation histogram range (log10 scale)';
            s.colortranslation.rangeMin = 1;
            doc.colortranslation.rangeMax = 'Maximum in histogram for computing color translation (log10 scale)';
            s.colortranslation.rangeMax = 5.5;
            doc.colortranslation.binIncrement = 'Resolution of histogram bins used for computing color translation (log10 scale)';
            s.colortranslation.binIncrement = 0.1;
            doc.colortranslation.minSamples = 'Minimum number of samples in a histogram bin for use in computing color translation';
            s.colortranslation.minSamples = 100;
            doc.colortranslation.channelMinimum = 'If set to [M1, M2, ...] trims channel i values below 10^Mi; otherwise drops those below 10^3';
            s.colortranslation.channelMinimum = {};
            doc.colortranslation.plot = 'Determines whether color translation plots should be created';
            s.colortranslation.plot = [];
            defaults('colortranslation.plot') = 'calibration.plot';
            doc.colortranslation.visiblePlots = 'If true, color translation plots are visible; otherwise, they are hidden for later saving';
            s.colortranslation.visiblePlots = [];
            defaults('colortranslation.visiblePlots') = 'calibration.visiblePlots';
            doc.colortranslation.plotPath = 'Location for color translation plots';
            s.colortranslation.plotPath = [];
            defaults('colortranslation.plotPath') = 'calibration.plotPath';
            doc.colortranslation.plotSize = 'Size (in inches) [X Y] for color translation figures';
            s.colortranslation.plotSize = [6 4];
            defaults('colortranslation.plotSize') = 'calibration.heatmapPlotSize';

            % Beads
            s.beads = struct(); doc.beads = struct();
            doc.beads.about = 'Settings controlling the interpretation of color calibration beads';
            doc.beads.catalogFileName = 'Location of bead catalog file';
            s.beads.catalogFileName = [fileparts(mfilename('fullpath')) '/../BeadCatalog.xlsx'];
            doc.beads.secondaryBeadChannel = 'For better distingishing low-a.u. ERF peaks: when set, segment ERF-channel peaks using the secondary channel instead of the ERF channel';
            s.beads.secondaryBeadChannel = [];
            doc.beads.peakThreshold = 'Manual minimum threshold for bead peaks; set automatically if empty';
            s.beads.peakThreshold = [];
            doc.beads.rangeMin = 'Minimum value considered for bead peaks (log scale: 10^rangeMin)';
            s.beads.rangeMin = 2;
            doc.beads.rangeMax = 'Maximum value considered for bead peaks (log scale: 10^rangeMax)';
            s.beads.rangeMax = 7;
            doc.beads.binIncrement = 'Resolution of histogram bins used for finding bead peaks (log10 scale)';
            s.beads.binIncrement = 0.02;
            doc.beads.beadModel = 'Model of beads that are being used. Should match an option in BeadCatalog.xlsx';
            s.beads.beadModel = 'SpheroTech RCP-30-5A';
            doc.beads.beadChannel = 'Laser/filter channel being used, defaults to FITC (MEFL). Should match an option in BeadCatalog.xlsx';
            s.beads.beadChannel = 'FITC';
            doc.beads.beadBatch = 'Batch of beads that are being used. If set, should match an option in BeadCatalog.xlsx';
            s.beads.beadBatch = [];
            doc.beads.forceFirstPeak = 'If set to N, lowest observed bead peak is forced to be interpreted as Nth peak';
            s.beads.forceFirstPeak = [];
            doc.beads.plot = 'When true, make diagnostic plots while computing bead unit calibration';
            s.beads.plot = [];
            defaults('beads.plot') = 'calibration.plot';
            doc.beads.visiblePlots = 'If true, bead unit calibration plots are visible; otherwise, they are hidden for later saving';
            s.beads.visiblePlots = [];
            defaults('beads.visiblePlots') = 'calibration.visiblePlots';
            doc.beads.plotPath = 'Location for bead unit calibration plots';
            s.beads.plotPath = [];
            defaults('beads.plotPath') = 'calibration.plotPath';
            doc.beads.plotSize = 'Size (in inches) [X Y] for bead unit calibration figures';
            s.beads.plotSize = [5 3.66];
            defaults('beads.plotSize') = 'calibration.graphPlotSize';
            doc.beads.validateAllChannels = 'If true, check all channels for likely bead problems; otherwise, check only ERF channel';
            s.beads.validateAllChannels = false;

            % Size Beads (for forward scatter calibration)
            s.sizebeads = struct(); doc.sizebeads = struct();
            doc.sizebeads.about = 'Settings controlling the interpretation of size calibration beads';
            doc.sizebeads.peakThreshold = 'Manual minimum threshold for size bead peaks; set automatically if empty';
            s.sizebeads.peakThreshold = [];
            doc.sizebeads.rangeMin = 'Minimum value considered for size bead peaks (log scale: 10^rangeMin)';
            s.sizebeads.rangeMin = 2;
            doc.beads.rangeMax = 'Maximum value considered for size bead peaks (log scale: 10^rangeMax)';
            s.sizebeads.rangeMax = 7;
            doc.sizebeads.binIncrement = 'Resolution of histogram bins used for finding size bead peaks (log10 scale)';
            s.sizebeads.binIncrement = 0.02;
            doc.sizebeads.beadModel = 'Model of size beads that are being used. Should match an option in BeadCatalog.xlsx';
            s.sizebeads.beadModel = 'SpheroTech PPS-6K';
            doc.sizebeads.beadChannel = 'Laser/filter channel being used, defaults to FSC. Should match an option in BeadCatalog.xlsx';
            s.sizebeads.beadChannel = 'FSC';
            doc.sizebeads.beadBatch = 'Batch of size beads that are being used. If set, should match an option in BeadCatalog.xlsx';
            s.sizebeads.beadBatch = [];
            doc.sizebeads.forceFirstPeak = 'If set to N, lowest observed size bead peak is forced to be interpreted as Nth peak';
            s.sizebeads.forceFirstPeak = [];
            doc.sizebeads.plot = 'When true, make diagnostic plots while computing size bead unit calibration';
            s.sizebeads.plot = [];
            defaults('sizebeads.plot') = 'calibration.plot';
            doc.sizebeads.visiblePlots = 'If true, size bead unit calibration plots are visible; otherwise, they are hidden for later saving';
            s.sizebeads.visiblePlots = [];
            defaults('sizebeads.visiblePlots') = 'calibration.visiblePlots';
            doc.sizebeads.plotPath = 'Location for size bead unit calibration plots';
            s.sizebeads.plotPath = [];
            defaults('sizebeads.plotPath') = 'calibration.plotPath';
            doc.sizebeads.plotSize = 'Size (in inches) [X Y] for size bead unit calibration figures';
            s.sizebeads.plotSize = [5 3.66];
            defaults('sizebeads.plotSize') = 'calibration.graphPlotSize';

            % OutputSettings migration
            doc.OutputSettings = struct();
            doc.OutputSettings.about = 'Settings controlling batch plotting';
            s.OutputSettings = struct();
            doc.OutputSettings.StemName = '';
            s.OutputSettings.StemName='';
            doc.OutputSettings.DeviceName = '';
            s.OutputSettings.DeviceName='';
            doc.OutputSettings.Description = '';
            s.OutputSettings.Description='';
            
            doc.OutputSettings.FixedInducerAxis = 'Set to fix limit [min max] of inducer count plot axis';
            s.OutputSettings.FixedInducerAxis = [];      % fixed -> [min max]
            doc.OutputSettings.FixedHistogramAxis = 'Set to fix limit [min max] of histogram count plot axis';
            s.OutputSettings.FixedHistogramAxis = [];
            doc.OutputSettings.FixedBinningAxis = 'Set to fix limit [min max] of binning variable plot axis';
            s.OutputSettings.FixedBinningAxis = [];
            doc.OutputSettings.FixedInputAxis = 'Set to fix limit [min max] of input plot axis';
            s.OutputSettings.FixedInputAxis =   [];
            doc.OutputSettings.FixedNormalizedInputAxis = 'Set to fix limit [min max] of normalized input plot axis';
            s.OutputSettings.FixedNormalizedInputAxis =   [];
            doc.OutputSettings.FixedOutputAxis = 'Set to fix limit [min max] of output plot axis';
            s.OutputSettings.FixedOutputAxis =  [];
            doc.OutputSettings.FixedNormalizedOutputAxis = 'Set to fix limit [min max] of normalized output plot axis';
            s.OutputSettings.FixedNormalizedOutputAxis =  [];
            doc.OutputSettings.FixedRatioAxis = 'Set to fix limit [min max] of ratio plot axis';
            s.OutputSettings.FixedRatioAxis =   [];
            doc.OutputSettings.FixedSNRAxis = 'Set to fix limit [min max] of signal-to-noise ratio plot axis';
            s.OutputSettings.FixedSNRAxis =   [];
            doc.OutputSettings.FixedDeltaSNRAxis = 'Set to fix limit [min max] of delta signal-to-noise ratio plot axis';
            s.OutputSettings.FixedDeltaSNRAxis =   [];
            doc.OutputSettings.ColorPlots = 'If true, color plots created';
            s.OutputSettings.ColorPlots = true;
            doc.OutputSettings.PlotPopulation = 'If true, population plots created';
            s.OutputSettings.PlotPopulation = true;
            doc.OutputSettings.PlotNormalized = 'If true, normalized plots created';
            s.OutputSettings.PlotNormalized = true;
            doc.OutputSettings.PlotNonnormalized = 'If true, nonnormalized plots created';
            s.OutputSettings.PlotNonnormalized = true;
            doc.OutputSettings.PlotEveryN = 'If true, plots every N';
            s.OutputSettings.PlotEveryN = 1;
            doc.OutputSettings.PlotTickMarks = 'If true, displays tick marks';
            s.OutputSettings.PlotTickMarks = false;
            doc.OutputSettings.FigureSize = 'Size (in inches) [X Y] for figures';
            s.OutputSettings.FigureSize = [5 3.66];
            doc.OutputSettings.csvfile = 'May be either a fid (file identifier) or a string';
            s.OutputSettings.csvfile = []; 
            
            % Plusminus preferences
            s.plusminus = struct(); 
            doc.plusminus = struct();
            doc.plusminus.about = 'Settings controlling plusminus plotting preferences';
            doc.plusminus.plotError = 'If true, plots error envelopes in plusminus comparison graphs';
            s.plusminus.plotError = false;
            
            % Histogram preferences
            s.histogram = struct(); 
            doc.histogram = struct();
            doc.histogram.about = 'Settings controlling histogram plotting preferences';
            doc.histogram.displayLegend = 'If true, displays legend in bin statistics graphs';
            s.histogram.displayLegend = true;
            
            % Excel wrapper preferences
            s.template = struct();
            doc.template = struct();
            doc.template.about = 'Settings controlling excel wrapper preferences';
            doc.template.displayErrors = 'If true, will display ALL of the TASBE warnings and errors from TemplateExtraction';
            s.template.displayErrors = false;
            
            % Last of all, bundle it in a cell array to return
            state = {s defaults doc};
        end
        
        function [out, default, doc] = setget(key,value,force)
            persistent state_stack; % contents: {{checkpoint {settings defaults documentation}}, ...}
            % if empty or reset, initialize
            is_reset = (nargin>0 && strcmp(key,'.reset'));
            if isempty(state_stack) || is_reset
                state_stack = {{'init' TASBEConfig.init()}}; 
                if is_reset, return; end;
            end;
            
            % unpack the current state stack
            settings = state_stack{end}{2}{1};
            defaults = state_stack{end}{2}{2};
            documentation = state_stack{end}{2}{3};
            
            % If there is no arguments, just return the current state for inspection
            if nargin==0, out = settings; default = defaults; doc = documentation; return; end
            % If the key is special, interpret it:
            % if the key is for checkpointing, then push/pop appropriately
            if strcmp(key,'.checkpoint'),
                % truncate to last named checkpoint, if any
                for i=1:numel(state_stack),
                    if(strcmp(state_stack{i}{1},value))
                        state_stack = {state_stack{1:(i-1)}}; % truncate
                        break;
                    end
                end
                % push a copy of the last onto the stack
                if ~isempty(state_stack)
                    state_stack{end+1} = {value state_stack{end}{2}};
                else
                    state_stack = {{'init' TASBEConfig.init()}}; 
                end
                out = state_stack{end}{2}{1};
                return
            end
            if strcmp(key,'.checkpoint_list'),
                out = cell(numel(state_stack),1);
                for i=1:numel(state_stack)
                    out{end-i+1} = state_stack{i}{1};
                end
                return
            end
            
            % Otherwise, go on to setting or getting
            if nargin<3, force = false; end
            
            keyseq = regexp(key, '\.', 'split');
            
            % nested access
            nest = cell(size(keyseq)); docNest = nest; 
            nest{1} = settings; docNest{1} = documentation;
            for i=1:(numel(keyseq)-1)
                if ~ismember(keyseq{i},fieldnames(nest{i})) || ~isstruct(nest{i}.(keyseq{i}))
                    nest{i}.(keyseq{i}) = struct();
                end
                if ~ismember(keyseq{i},fieldnames(docNest{i})) || ~isstruct(docNest{i}.(keyseq{i}))
                    docNest{i}.(keyseq{i}) = struct();
                end
                nest{i+1} = nest{i}.(keyseq{i});
                docNest{i+1} = docNest{i}.(keyseq{i});
            end
            
            if ~isempty(value) || force
                % warn if the key is not already in the field names
                if ~ismember(keyseq{end},fieldnames(nest{end}))
                    TASBESession.warn('TASBEConfig','UnknownSetting','Setting previously unknown configuration "%s"',key);
                end
                % set the value
                nest{end}.(keyseq{end}) = value;
                % propagate up the chain
                for i=1:(numel(keyseq)-1)
                    nest{end-i}.(keyseq{end-i}) = nest{end-i+1};
                end
                settings = nest{1};
                state_stack{end}{2}{1} = settings; % remember it for next time
            end
            
            % Finally, set outputs
            if ismember(keyseq{end},fieldnames(nest{end}))
                out = nest{end}.(keyseq{end});
            else
                out = [];
            end
            if ismember(keyseq{end},fieldnames(docNest{end}))
                doc = docNest{end}.(keyseq{end});
            else
                doc = 'No documentation available';
            end
            
            if(defaults.isKey(key)), default = defaults(key); else default = 'no default'; end;
        end
    end
    
    methods(Static)
        % Set a value for this key
        function out = set(key,value)
            out = TASBEConfig.setget(key,value);
        end
        
        % Get a value for this key, not following its default chain. If no value is set, try to use the 2nd argument
        function out = getexact(key,default)
            % try a get
            out = TASBEConfig.setget(key,[]);
            if isempty(out)
                % if empty, try to set to default
                if nargin>=2
                    out = TASBEConfig.setget(key,default);
                else
                    error('TASBEConfig', 'NoDefault', 'Requested non-existant setting without default: %s',key);
                end
            end
        end
        
        % Get the first defined value in a sequence of possibilities
        function out = getseq(varargin)
            for i=1:numel(varargin)
                try
                    out = TASBEConfig.get(varargin{i});
                    return
                catch e % ignore error and continue
                end
            end
            error('TASBEConfig', 'NoPreference', 'Couldn''t get any preference in sequence: %s',[varargin{:}]);
        end
        
        % Get a value for this key, possibly via default
        function out = get(key)
            persistent defaults
            if isempty(defaults), [~, defaults] = TASBEConfig.list(); end;
            
            current = key;
            while current
                try
                    out = TASBEConfig.getexact(current);
                    return
                catch e % on miss, try to follow preference path
                    try
                        current = defaults(current);
                    catch e
                        error('TASBEConfig', 'NoPreference', 'Couldn''t get any preference for: %s',key);
                    end
                end
            end
        end
        
        % Check if there is a value set for this key or one of its defaults
        function TF = isSet(key)
            try
                TASBEConfig.get(key);
                TF = true;
            catch % error means not set
                TF = false;
            end
        end
        
        % Remove any value set for key
        function clear(key)
            TASBEConfig.setget(key,[],true);
        end
        
        % Reset the entire TASBEConfig state
        function reset()
            TASBEConfig.setget('.reset');
        end
        
        % Return a list of all of the checkpoint names
        function [top, chain] = checkpoints() 
            chain = TASBEConfig.setget('.checkpoint_list');
            top = chain{1};
        end
        
        % 
        function old = checkpoint(name) 
            old = TASBEConfig.setget('.checkpoint',name);
        end
        
        % Get help on a particular key
        function text = help(key)
            if nargin==0
                [val,def,doc] = TASBEConfig.list();
                key = '';
            else
                [val,def,doc] = TASBEConfig.setget(key,[]);
            end
            
            if(isstruct(val))
                fieldnameset = fieldnames(val);
                keydoc = '';
                maxname = max(cellfun(@numel,fieldnameset));
                for i = 1:numel(fieldnameset),
                    keydoc = [keydoc sprintf('\n  %s',fieldnameset{i})];
                    if(isstruct(val.(fieldnameset{i}))), 
                        keydoc = [keydoc sprintf('\t\t[family]')]; 
                    else
                        try 
                            about = doc.(fieldnameset{i}); 
                            spacer = repmat(' ',1, 4 + maxname - numel(fieldnameset{i}));
                            keydoc = [keydoc sprintf('%s%s',spacer,about)];
                        catch e, 
                            % leave it as was, without documentation
                        end;
                    end;
                end
                try about = doc.about; catch e, about = 'No documentation available'; end;
                text = sprintf('Configuration family: %s\n%s\nKeys: %s',key,about,keydoc);
            else
                text = sprintf('Configuration: %s\n%s\nDefaults to: %s\nCurrent value: %s\n',key,doc,def,val);
            end
        end
        
        % Get a list of all keys
        function [settings, defaults, documentation] = list()
            [settings, defaults, documentation] = TASBEConfig.setget();
        end
        
        % Transform all non-default settings into a JSON object
        function string = to_json()
            string = savejson('',TASBEConfig.list());
        end
        
        function load_from_json(string)
            json_object = loadjson(string);
            TASBEConfig.set_config_from_object('', json_object);
        end
        
        function set_config_from_object(prefix, struct)
            fields = fieldnames(struct);
            for i=1:numel(fields);
                value = struct.(fields{i});
                if(isstruct(value)), % set with substructure
                    TASBEConfig.set_config_from_object([prefix fields{i} '.'], value);
                else % set this value
                    TASBEConfig.set([prefix fields{i}], value);
                end
            end
        end
    end
end
