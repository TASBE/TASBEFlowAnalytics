% get_channel_results is a getter function for the ExperimentResults class
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.


function [means stds stdofmeans stdofstds] = get_channel_results(ER,name)
        channel = getChannel(ER.AnalysisParameters, name);
        which = indexof({ER.Means{:,1}},channel);
        means = ER.Means{which,2};
        which = indexof({ER.StandardDevs{:,1}},channel);
        stds = ER.StandardDevs{which,2};
        which = indexof({ER.StdOfMeans{:,1}},channel);
        stdofmeans = ER.StdOfMeans{which,2};
        which = indexof({ER.StdOfStandardDevs{:,1}},channel);
        stdofstds = ER.StdOfStandardDevs{which,2};
