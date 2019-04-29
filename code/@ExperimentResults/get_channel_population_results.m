% get_channel_population_results is a getter function for the ExperimentResults
% class, obtains the bulk statistics of the expression population, including histograms
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
function [means stds histograms excluded stdofmeans stdofstds stdofhistograms stdofexcluded] = get_channel_population_results(ER,name)
% function [means stds histograms excluded stdofmeans stdofstds stdofhistograms stdofexcluded] = get_channel_population_results(ER,name)

        channel = getChannel(ER.AnalysisParameters, name);
        which = indexof({ER.PopMeans{:,1}},channel);
        means = ER.PopMeans{which,2};
        which = indexof({ER.PopStandardDevs{:,1}},channel);
        stds = ER.PopStandardDevs{which,2};
        which = indexof({ER.Histograms{:,1}},channel);
        histograms = ER.Histograms{which,2};
        which = indexof({ER.Excluded{:,1}},channel);
        excluded = ER.Excluded{which,2};
        
        which = indexof({ER.StdOfPopMeans{:,1}},channel);
        stdofmeans = ER.StdOfPopMeans{which,2};
        which = indexof({ER.StdOfPopStandardDevs{:,1}},channel);
        stdofstds = ER.StdOfPopStandardDevs{which,2};
        which = indexof({ER.StdOfHistograms{:,1}},channel);
        stdofhistograms = ER.StdOfHistograms{which,2};
        which = indexof({ER.StdOfExcluded{:,1}},channel);
        stdofexcluded = ER.StdOfExcluded{which,2};
