% Constructor of the UnitTranslation class with the following properties:
% It is assumed that an array is used to associated this with channels
%    beadtype    % string, right now always RCP-30-5A
%    k_ERF      % ERF = k_ERF*ERF_channel_AU
%    first_peak  % what is the first RCP-30-5A peak visible?
%    fit_error   % residual from the linear fit
%    peak_sets   % all peaks from all channels, for logging and comparison
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
    
function UT = UnitTranslation(beadtype, k_ERF, first_peak, fit_error, peak_sets)
  if nargin == 0
      UT.beadtype = '';
      UT.k_ERF = 0;
      UT.first_peak = 0;
      UT.fit_error = 0;
      UT.peak_sets = {};
  elseif nargin == 5
      UT.beadtype = beadtype;
      UT.k_ERF = k_ERF;
      UT.first_peak = first_peak;
      UT.fit_error = fit_error;
      UT.peak_sets = peak_sets;
  end
  UT=class(UT,'UnitTranslation');

