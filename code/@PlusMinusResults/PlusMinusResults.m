% Constructor of PlusMinusResults class with the following properties:
% PlusResults             % results for plus experiment, including experiment and analysis parameters
% MinusResults            % results for plus experiment, including experiment and analysis parameters
% MeanRatio               % mean ratio of all valid ratios
% BinCounts               % table of n_bins
% Valid                   % table of n_bins x n_inductions
% Ratios                  % table of n_bins x n_inductions
% InMeans                 % table of n_bins X n_inductions x {plus,minus}
% OutMeans                % table of n_bins X n_inductions x {plus,minus}
% OutStandardDevs         % table of n_bins X n_inductions x {plus,minus}
% OutStdOfMeans           % table of n_bins X n_inductions x {plus,minus}
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
       
function PMR = PlusMinusResults(PlusResults, MinusResults, MeanRatio, BinCounts, Valid, Ratios, InMeans, InStandardDevs, InStdOfMeans, OutMeans, OutStandardDevs, OutStdOfMeans)
        PMR.PlusResults = ExperimentResults();
        PMR.MinusResults = ExperimentResults();
        PMR.MeanRatio = [];
        PMR.BinCounts =[];
        PMR.Valid = [];
        PMR.Ratios = [];
        PMR.InMeans = [];
        PMR.InStandardDevs = [];
        PMR.InStdOfMeans = [];
        PMR.OutMeans = [];
        PMR.OutStandardDevs = [];
        PMR.OutStdOfMeans = [];
   if nargin > 0
        PMR.PlusResults = PlusResults;
        PMR.MinusResults = MinusResults;
        PMR.MeanRatio = MeanRatio;
        PMR.BinCounts =BinCounts;
        PMR.Valid = Valid;
        PMR.Ratios =Ratios;
        PMR.InMeans = InMeans;
        PMR.InStandardDevs = InStandardDevs;
        PMR.InStdOfMeans = InStdOfMeans;
        PMR.OutMeans = OutMeans;
        PMR.OutStandardDevs = OutStandardDevs;
        PMR.OutStdOfMeans = OutStdOfMeans;
    end;
