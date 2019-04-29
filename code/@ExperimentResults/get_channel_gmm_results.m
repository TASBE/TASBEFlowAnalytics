% get_channel_gmm_results is a getter function for the ExperimentResults
% class, obtains the gaussian mixture model (population 'component') stats
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [means stds weights stdofmeans stdofstds stdofweights] = get_channel_gmm_results(ER,name)
% function [means stds weights stdofmeans stdofstds stdofweights] = get_channel_gmm_results(ER,name)
        channel = getChannel(ER.AnalysisParameters, name);
        which = indexof({ER.PopComponentMeans{:,1}},channel);
        means = ER.PopComponentMeans{which,2};
        which = indexof({ER.PopComponentStandardDevs{:,1}},channel);
        stds = ER.PopComponentStandardDevs{which,2};
        which = indexof({ER.PopComponentWeights{:,1}},channel);
        weights = ER.PopComponentWeights{which,2};
        
        which = indexof({ER.StdOfPopComponentMeans{:,1}},channel);
        stdofmeans = ER.StdOfPopComponentMeans{which,2};
        which = indexof({ER.StdOfPopComponentStandardDevs{:,1}},channel);
        stdofstds = ER.StdOfPopComponentStandardDevs{which,2};
        which = indexof({ER.StdOfPopComponentWeights{:,1}},channel);
        stdofweights = ER.StdOfPopComponentWeights{which,2};
