% Constructor for the SampleResults class with the following properties:
% Experiment              % handle to object, which includes induction conditions
% File                    % name of the file this came from
% AnalysisParameters      % handle to object which was used for producing these results
% Constitutive-correlated statistics
% BinCounts               % table of n_bins
% Means                   % table of n_bins X n_channels
% StandardDevs            % table of n_bins X n_channels
% FractionActive          % array of n_bins
% PlasmidEstimates        % array of n_bins
% PlasmidModel            % PlasmidExpressionModel for correlated
% Bulk statistics
% Histograms              % table of n_bins x n_channels Note: Histograms(:,constitutive) = BinCounts
% PopMeans                % array of n_channels
% PopStandardDevs         % array of n_channels
% PopPeaks                % cell array of n_channels, each containing a sorted list of histogram peaks
% Excluded                % array of n_channels: number of data points in sample not within bin range for each channel
% NonExpressing           % events excluded on low end from all channels
% k-component gaussian model
% PopComponentMeans       % table of n_components X n_channels
% PopComponentStandardDevs % table of n_components X n_channels
% PopComponentWeights     % table of n_components X n_channels
% Frac of points in on group (determined by threshold)
% Frac of points in off group (determined by threshold)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function SR = SampleResults(Experiment, File, AnalysisParameters, BinCounts, Means, StandardDevs, ...
    FractionActive, PlasmidEstimates, PlasmidModel, ...
    Histograms, PopMeans, PopStandardDevs, PopPeaks, Excluded, NonExpressing, ...
    PopComponentMeans, PopComponentStandardDevs, PopComponentWeights, on_frac, off_frac)
        SR.version = tasbe_version();
        SR.Experiment = Experiment();
        SR.File = File;
        SR.AnalysisParameters = AnalysisParameters();
        % Constitutive-correlated statistics
        SR.BinCounts =[];
        SR.Means = [];
        SR.StandardDevs = [];
        SR.FractionActive = [];
        SR.PlasmidEstimates = [];
        SR.PlasmidModel = PlasmidExpressionModel();    
        % Bulk statistics
        SR.Histograms = [];
        SR.PopMeans = [];
        SR.PopStandardDevs = [];
        SR.PopPeaks = [];
        SR.Excluded = [];
        SR.NonExpressing = [];
        % k-component gaussian model
        SR.PopComponentMeans = [];
        SR.PopComponentStandardDevs = [];
        SR.PopComponentWeights = [];
        SR.on_frac = 1;
        SR.off_frac = 0;
   if nargin > 0
        SR.Experiment = Experiment;
        SR.File = File;
        SR.AnalysisParameters = AnalysisParameters;
        % Constitutive-correlated statistics
        SR.BinCounts = BinCounts;
        SR.Means = Means;
        SR.StandardDevs = StandardDevs;
        SR.FractionActive = FractionActive;
        SR.PlasmidEstimates = PlasmidEstimates;
        SR.PlasmidModel = PlasmidModel;
        % Bulk statistics
        SR.Histograms = Histograms;
        SR.PopMeans = PopMeans;
        SR.PopStandardDevs = PopStandardDevs;
        SR.PopPeaks = PopPeaks;
        SR.Excluded = Excluded;
        SR.NonExpressing = NonExpressing;
        % k-component gaussian model
        SR.PopComponentMeans = PopComponentMeans;
        SR.PopComponentStandardDevs = PopComponentStandardDevs;
        SR.PopComponentWeights = PopComponentWeights;
        if nargin > 18
            SR.on_frac = on_frac;
            SR.off_frac = off_frac;
        end
   end
    
  %FY: I intentionally commented out the following class statement. 
  %There are no methods for this class thus it is as good as a struct.
  
  %SR=class(SR,'SampleResults');  

