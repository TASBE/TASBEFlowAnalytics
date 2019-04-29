% EXPERIMENTRESULTS summarizes the results of the experiment across all
% induction conditions and replications, contains the following properties:
%   
%         Experiment              % handle to object, which includes induction conditions
%         AnalysisParameters      % handle to object which was used for producing these results
%         % Means across all replications
%         ReplicateCounts         % table of n_bins x n_inductions
%         BinCounts               % table of n_bins x n_inductions
%         Means                   % table of {Channel,n_bins X n_inductions x n_channels} mean over replications
%         StandardDevs            % table of {Channel,n_bins X n_inductions; ...} mean over replications
%         PlasmidEstimates        % array of n_bins x n_inductions
%         FractionActive          % array of n_bins x n_inductions
%         PopMeans                % table of {Channel,n_inductions}
%         PopStandardDevs         % table of {Channel,n_inductions}
%         Histograms              % table of {Channel,n_bins x n_inductions}  (normalized to sample size)
%         Excluded                % table of {Channel,n_inductions}           (normalized to sample size)
%         NonExpressing           % table of n_inductions                     (normalized to sample size)
%         PopComponentMeans       % table of {Channel,n_components x n_inductions}
%         PopComponentStandardDevs % table of {Channel,n_components x n_inductions}
%         PopComponentWeights     % table of {Channel,n_components x n_inductions}
%
%         % Std devs across all replications
%         StdOfMeans              % table of {Channel,n_bins X n_inductions; ...} std over replications
%         StdOfStandardDevs       % table of {Channel,n_bins X n_inductions; ...} std over replications
%         StdOfPlasmidEstimates   % array of n_bins x n_inductions
%         StdOfFractionActive     % array of n_bins x n_inductions
%         StdOfPopMeans           % table of {Channel,n_inductions}
%         StdOfPopStandardDevs    % table of {Channel,n_inductions}
%         StdOfHistograms         % table of {Channel,n_bins x n_inductions}  (normalized to sample size)
%         StdOfExcluded           % table of {Channel,n_inductions}           (normalized to sample size)
%         StdOfNonExpressing      % table of n_inductions                     (normalized to sample size)
%         StdOfPopComponentMeans      % table of {Channel,n_components x n_inductions}
%         StdOfPopComponentStandardDevs % table of {Channel,n_components x n_inductions}
%         StdOfPopComponentWeights    % table of {Channel,n_components x n_inductions}
%         on_fracMean                 % arithmetic mean of frac of points in on group between replicates (determined by threshold)
%         off_fracMean                % arithmetic mean of frac of points in off group between replicates (determined by threshold)
%         on_fracStd                 % arithmetic std of frac of points in on group between replicates (determined by threshold)
%         off_fracStd                % arithmetic std of frac of points in off group between replicates (determined by threshold)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function ER = ExperimentResults(Exp, AP, ReplicateCounts, BinCounts, ...
                        Means, StandardDevs, PlasmidEstimates, FractionActive, ...
                        PopMeans, PopStandardDevs, Histograms, Excluded, NonExpressing, ...
                        PopComponentMeans, PopComponentStandardDevs, PopComponentWeights, ...
                        StdOfMeans, StdOfStandardDevs, StdOfPlasmidEstimates, StdOfFractionActive, ...
                        StdOfPopMeans, StdOfPopStandardDevs, StdOfHistograms, StdOfExcluded, StdOfNonExpressing, ...
                        StdOfPopComponentMeans, StdOfPopComponentStandardDevs, StdOfPopComponentWeights, on_fracMean, off_fracMean, on_fracStd, off_fracStd)
                ER.version = tasbe_version();
            ER.Experiment = Experiment();
                ER.AnalysisParameters = AnalysisParameters();
                ER.ReplicateCounts = [];
                ER.BinCounts = [];
                channels = cell(1, 2);
                channels{1,1} = Channel();
                channels{1,2} = [];
                ER.Means = channels;
                ER.StandardDevs = channels;
                ER.PlasmidEstimates = [];
                ER.FractionActive = [];
                ER.PopMeans = channels;
                ER.PopStandardDevs = channels;
                ER.Histograms = channels;
                ER.Excluded = channels;
                ER.NonExpressing = [];
                ER.PopComponentMeans = channels;
                ER.PopComponentStandardDevs = channels;
                ER.PopComponentWeights = channels;
                ER.StdOfMeans = channels;
                ER.StdOfStandardDevs = channels;
                ER.StdOfPlasmidEstimates = [];
                ER.StdOfFractionActive = [];
                ER.StdOfPopMeans = channels;
                ER.StdOfPopStandardDevs = channels;
                ER.StdOfHistograms = channels;
                ER.StdOfExcluded = channels;
                ER.StdOfNonExpressing = [];
                ER.StdOfPopComponentMeans = channels;
                ER.StdOfPopComponentStandardDevs = channels;
                ER.StdOfPopComponentWeights = channels;
            if nargin > 0
                ER.Experiment = Exp;
                ER.AnalysisParameters = AP;
                ER.ReplicateCounts = ReplicateCounts;
                ER.BinCounts = BinCounts;
                ER.Means = Means;
                ER.StandardDevs = StandardDevs;
                ER.PopMeans = PopMeans;
                ER.PopStandardDevs = PopStandardDevs;
                ER.Histograms = Histograms;
                ER.Excluded = Excluded;
                ER.NonExpressing = NonExpressing;
                ER.PopComponentMeans = PopComponentMeans;
                ER.PopComponentStandardDevs = PopComponentStandardDevs;
                ER.PopComponentWeights = PopComponentWeights;
                ER.PlasmidEstimates = PlasmidEstimates;
                ER.FractionActive = FractionActive;
                ER.StdOfMeans = StdOfMeans;
                ER.StdOfStandardDevs = StdOfStandardDevs;
                ER.StdOfPlasmidEstimates = StdOfPlasmidEstimates;
                ER.StdOfFractionActive = StdOfFractionActive;
                ER.StdOfPopMeans = StdOfPopMeans;
                ER.StdOfPopStandardDevs = StdOfPopStandardDevs;
                ER.StdOfHistograms = StdOfHistograms;
                ER.StdOfExcluded = StdOfExcluded;
                ER.StdOfNonExpressing = [];
                ER.StdOfPopComponentMeans = StdOfPopComponentMeans;
                ER.StdOfPopComponentStandardDevs = StdOfPopComponentStandardDevs;
                ER.StdOfPopComponentWeights = StdOfPopComponentWeights;
                if nargin > 28
                    ER.on_fracMean = on_fracMean;
                    ER.off_fracMean = off_fracMean;
                    ER.on_fracStd = on_fracStd;
                    ER.off_fracStd = off_fracStd;
                end
            end
            
            ER=class(ER,'ExperimentResults');
