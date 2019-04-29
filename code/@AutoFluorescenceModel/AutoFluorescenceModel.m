% AutoFluorescenceModel constructs an AutoFluorescenceModel object for its
% channel using statistics from the provided data.
%
% It is assumed that an array is used to associate this with channels
%  af_mean          % linear mean in arbitrary FACS units
%  af_std           % linear std.dev. in arbitrary FACS units
%  af_mean_ERF=[] v % linear mean in ERF, set when ERFized
%  af_std_ERF=[]    % linear std.dev. in ERF, set when ERFized
%  n                % number of points used in this computation
%  channel          % pointer to the channel for which this is a model
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function AFM = AutoFluorescenceModel(channel,data)
    if nargin == 0
        AFM.channel = [];
        AFM.af_mean = 0;
        AFM.af_std = 0;
        AFM.n = 0;
    elseif nargin >= 2
        AFM.channel = channel;
        % to exclude outliers, drop top and bottom fractions of data
        sorted = sort(data);
        dropFraction = TASBEConfig.get('autofluorescence.dropFraction');
        dropsize = ceil(numel(sorted)*dropFraction);
        trimmed = sorted(dropsize:(numel(sorted)-dropsize));
        % compute statistics
        AFM.af_mean = mean(trimmed);
        AFM.af_std = std(trimmed);
        AFM.n = numel(trimmed);
    else
        TASBESession.error('AutoFluorescence','MissingArgument','Autofluorescence Model constructor requires two arguments');
    end
    % make placeholders for eventual ERF conversion information
    AFM.af_mean_ERF = [];
    AFM.af_std_ERF = [];
    AFM=class(AFM,'AutoFluorescenceModel');
       
