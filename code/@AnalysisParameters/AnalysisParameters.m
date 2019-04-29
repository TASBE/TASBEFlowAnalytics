% ANALYSISPARAMETERS contains parameters that are necessary for analysis of experiment
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function AP = AnalysisParameters(bins, ChannelLabels, ERFPerPlasmid)
   AP.version = tasbe_version();
   AP.ERFPerPlasmid = 1000;                % default value 1000
   AP.use_autofluorescence = false;         % if false, try to remove autofluorescence from results
   AP.min_valid_count = 100;                % minimum number of cells per bin to be allowed
   AP.min_fraction_active = 0.9;            % minimum fraction active to be considered valid in a bin
   AP.pem_drop_threshold = 5;               % below what level should the plasmid model ignore data?
   AP.num_gaussian_components = 2;           % How many major populations of cells are expected? (default 1 active, 1 inactive)

   if nargin == 0 
       AP.bins = BinSequence();                          % a BinSequence object
       AP.ChannelLabels = {'',Channel()};        % e.g. {'constitutive', bluechannel; 'input', redchannel}       
   elseif nargin >= 2
       AP.bins = bins;                          % a BinSequence object
       AP.ChannelLabels = ChannelLabels;        % e.g. {'constitutive', bluechannel; 'input', redchannel}
   end
   
   if nargin > 2, AP.ERFPerPlasmid = ERFPerPlasmid; end
   AP=class(AP,'AnalysisParameters');


