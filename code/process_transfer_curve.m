% PROCESS_TRANSFER_CURVE analyzes the data for transfer curve analysis.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [results, sampleresults] = process_transfer_curve( colorModel, experiment, analysisParams)

data = read_data( colorModel, experiment, analysisParams);

sampleresults = process_data(colorModel,experiment,analysisParams, data);

results = summarize_data(colorModel,experiment,analysisParams,sampleresults);
