% Copyright (C) 2011 - 2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE Flow Analytics distribution's top directory.
%
% This file is part of the TASBE Flow Analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the BBN Flow Cytometry
% package distribution's top directory.

classdef TASBEConfig
    methods(Static,Hidden)
        function [s, defaults] = init()
            s = struct();
            defaults = containers.Map();
            
            % Testing variables
            s.testing.fakeFigureSaves = 0;
            
            % Matlab GMdistribution
            
            % Generic flow data analysis
            s.flow.rangeMin = 0;                           % bin minimum (log10 scale)
            s.flow.rangeMax = 7;                           % bin maximum (log10 scale)
            s.flow.outputPointCloud = false;               % if true, output pointcloud for each calibrated read
            s.flow.pointCloudPath = 'CSV/';                % location for pointcloud outputs

            % generic plots
            s.plotPath = 'plots/';          % where should any plot go?
            s.visiblePlots = false;         % Should plots be visible, or just created?
            s.graphPlotSize = [6 4];        % What size (in inches) should data graph figures be?
            s.heatmapPlotSize = [5 5];      % What size (in inches) should scatter/heatmap figures be?
            s.largeOutliers = false;        % Should outliers in heatmap figures be large, for output in small figures?
            s.heatmapPlotType = 'image';    % Should contour plots be used instead of heatmap images?
            
            % supporting plots, i.e., those supporting the transformation of raw data into processed data, like autofluorescence, compensation, units
            s.supporting = struct();
            s.supporting.plot = true;           % make plots as a side effect of computing color models, etc.
            s.supporting.visiblePlots = [];     % should supporting plots be visible, or just created?
            defaults('supporting.visiblePlots') = 'visiblePlots';
            s.supporting.plotPath = [];         % where should supporting plots go?
            defaults('supporting.plotPath') = 'plotPath';
            s.supporting.graphPlotSize = [];    % What size (in inches) should supporting data graph figures be?
            defaults('supporting.graphPlotSize') = 'graphPlotSize';
            s.supporting.heatmapPlotSize = [];  % What size (in inches) should supporting scatter/heatmap figures be?
            defaults('supporting.heatmapPlotSize') = 'heatmapPlotSize';
            
            % Gating
%             s.gating = struct();
%             s.gating.fractionFromExtrema = 0.95;    % Fraction of range considered not-saturated and thus included in gating
%             s.gating.saturationWarning = 0.3;       % Warn about saturation distorting gates if fraction non-extrema less than this level
%             s.gating.numComponents = 2;             % number of gaussian components searched for
%             s.gating.rankedComponents = 1;          % Array of which components will be selected, in order of tightness
%             s.gating.deviations = 2.0;              % number of standard deviations out the gaussian that will be allowed
%             s.gating.tightening = [];               % If set, amount that selected components are further tightened (range: [0,1])
%             s.gating.plot = [];                     % Should a gating plot be created?
%             s.gating.showNonselected = true;        % Should plot show only the selected component(s), or all of them?
%             defaults('gating.plot') = 'supporting.plot';
%             s.gating.visiblePlots = [];             % should gating plot be visible, or just created?
%             defaults('gating.visiblePlots') = 'supporting.visiblePlots';
%             s.gating.plotPath = [];                 % where should gating plot go?
%             defaults('gating.plotPath') = 'supporting.plotPath';
%             s.gating.plotSize = [];                 % What size (in inches) should gating plot be?
%             defaults('gating.plotSize') = 'supporting.heatmapPlotSize';
            
            % Autofluorescence removal
%             s.autofluorescence = struct();
%             s.autofluorescence.dropFractions = 0.025;   % What fraction of the extrema should be dropped before computing autofluorescence?
%             s.autofluorescence.plot = [];               % Should an autofluorescence plot be created?
%             defaults('autofluorescence.plot') = 'supporting.plot';
%             s.autofluorescence.visiblePlots = [];       % should autofluorescence plot be visible, or just created?
%             defaults('autofluorescence.visiblePlots') = 'supporting.visiblePlots';
%             s.autofluorescence.plotPath = [];           % where should autofluorescence plot go?
%             defaults('autofluorescence.plotPath') = 'supporting.plotPath';
%             s.autofluorescence.plotSize = [];           % What size (in inches) should autofluorescence plot be?
%             defaults('autofluorescence.plotSize') = 'supporting.graphPlotSize';
            
            % Spectral bleed compensation
            s.compensation = struct();
            s.compensation.minimumDrivenLevel = 1e2;    % uniformly ignore all less than this level of a.u. 
            s.compensation.maximumDrivenLevel = Inf;     % uniformly ignore all greater than this level of a.u. 
            s.compensation.minimumBinCount = 10;        % ignore bins with less than this many elements
            s.compensation.highBleedWarning = 0.1;      % Warn about high bleed at this level
            s.compensation.plot = [];                   % Should compensation plots be created?
            defaults('compensation.plot') = 'supporting.plot';
            s.compensation.visiblePlots = [];           % should compensation plot be visible, or just created?
            defaults('compensation.visiblePlots') = 'supporting.visiblePlots';
            s.compensation.plotPath = [];               % where should compensation plot go?
            defaults('compensation.plotPath') = 'supporting.plotPath';
            s.compensation.plotSize = [];               % What size (in inches) should compensation figure be?
            defaults('compensation.plotSize') = 'supporting.heatmapPlotSize';
            
            % Beads
%             s.beads = struct();
%             s.beads.catalogFileName = 'BeadCatalog.xls';    % Where is the catalog file?
%             s.beads.peakThreshold = [];                     % Manual minimum threshold for peaks; set automatically if empty
%             s.beads.rangeMin = 2;                           % bin minimum (log10 scale)
%             s.beads.rangeMax = 7;                           % bin maximum (log10 scale)
%             s.beads.binIncrement = 0.02;                    % resolution of binning
%             s.beads.plot = [];                              % Should an autofluorescence plot be created?
%             s.beads.forceFirstPeak = [];                    % If set to N, lowest observed peak is forced to be batch to Nth peak
%             defaults('beads.plot') = 'supporting.plot';
%             s.beads.visiblePlots = [];                      % should autofluorescence plot be visible, or just created?
%             defaults('beads.visiblePlots') = 'supporting.visiblePlots';
%             s.beads.plotPath = [];                          % where should autofluorescence plot go?
%             defaults('beads.plotPath') = 'supporting.plotPath';
%             s.beads.plotSize = [];                          % What size (in inches) should autofluorescence plot be?
%             defaults('beads.plotSize') = 'supporting.graphPlotSize';
            
            % Color translation
%             s.colortranslation = struct();
%             s.colortranslation.rangeMin = 3;                % bin minimum (log10 scale), universal minimum trim
%             s.colortranslation.rangeMax = 5.5;              % bin maximum (log10 scale)
%             s.colortranslation.binIncrement = 0.1;          % resolution of binning
%             s.colortranslation.minSamples = 100;            % How many samples are needed for a bin's data to be used?
%             s.colortranslation.trimMinimum = {};            % If set, trims individual channels via {{Channel,log10(min)} ...}
%             s.colortranslation.plot = [];                   % Should an autofluorescence plot be created?
%             defaults('colortranslation.plot') = 'supporting.plot';
%             s.colortranslation.visiblePlots = [];           % should autofluorescence plot be visible, or just created?
%             defaults('colortranslation.visiblePlots') = 'supporting.visiblePlots';
%             s.colortranslation.plotPath = [];               % where should autofluorescence plot go?
%             defaults('colortranslation.plotPath') = 'supporting.plotPath';
%             s.colortranslation.plotSize = [];               % What size (in inches) should autofluorescence plot be?
%             defaults('colortranslation.plotSize') = 'supporting.heatmapPlotSize';
        end
        
        function out = setget(key,value,force)
            persistent settings;
            if isempty(settings), settings = TASBEConfig.init(); end;
            
            % If there is no arguments, just return the whole thing for inspection
            if nargin==0, out = settings; return; end
            % if the key is '.reset', then reset
            if strcmp(key,'.reset'), settings = TASBEConfig.init(); return; end; 
            if nargin<3, force = false; end
            
            keyseq = regexp(key, '\.', 'split');
            
            % nested access
            nest = cell(size(keyseq));
            nest{1} = settings;
            for i=1:(numel(keyseq)-1)
                if ~ismember(keyseq{i},fieldnames(nest{i})) || ~isstruct(nest{i}.(keyseq{i}))
                    nest{i}.(keyseq{i}) = struct();
                end
                nest{i+1} = nest{i}.(keyseq{i});
            end
            
            if ~isempty(value) || force
                nest{end}.(keyseq{end}) = value;
                % propagate up the chain
                for i=1:(numel(keyseq)-1)
                    nest{end-i}.(keyseq{end-i}) = nest{end-i+1};
                end
                settings = nest{1};
            end
            
            if ismember(keyseq{end},fieldnames(nest{end}))
                out = nest{end}.(keyseq{end});
            else
                out = [];
            end
        end
    end
    
    methods(Static)
        function out = set(key,value)
            out = TASBEConfig.setget(key,value);
        end
        
        function out = getexact(key,default)
            % try a get
            out = TASBEConfig.setget(key,[]);
            if isempty(out)
                % if empty, try to set to default
                if nargin>=2
                    out = TASBEConfig.setget(key,default);
                else
                    error('Requested non-existing setting without default: %s',key);
                end
            end
        end
        
        % Get the first defined in a sequence
        function out = getseq(varargin)
            for i=1:numel(varargin)
                try
                    out = TASBEConfig.get(varargin{i});
                    return
                catch e % ignore error and continue
                end
            end
            error('Couldn''t get any preference in sequence: %s',[varargin{:}]);
            
        end
        function out = get(key)
            persistent defaults
            if isempty(defaults), [ignored defaults] = TASBEConfig.init(); end;
            
            current = key;
            while current
                try
                    out = TASBEConfig.getexact(current);
                    return
                catch e % on miss, try to follow preference path
                    try
                        current = defaults(current);
                    catch e
                        error('Couldn''t get any preference for: %s',key);
                    end
                end
            end
        end
        
        function TF = isSet(key)
            try
                TASBEConfig.get(key);
                TF = true;
            catch % error means not set
                TF = false;
            end
        end
        function clear(key)
            TASBEConfig.setget(key,[],true);
        end
        function reset()
            TASBEConfig.setget('.reset');
        end
        
        function out = list()
            out = TASBEConfig.setget();
        end
    end
end
