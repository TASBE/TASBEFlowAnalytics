% Constructor of NoiseModel structure with the following properties:
% noisemin    % array of min noise/channel (which is typically the highest)
% noisemean   % array of std noise/channel
% noisestd    % array of noise std/channel (for checking validity)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function NM = NoiseModel(mins, means, stds, detailcounts, detailmeans, detailstds)
  NM.noisemin = mins;
  NM.noisemean = means;
  NM.noisestd = stds;
  NM.detailcounts = detailcounts;
  NM.detailmeans = detailmeans;
  NM.detailstds = detailstds;
  % this line is commented intentionally. If we ever to decide
  % methods for this class then we should turn this into a real
  % class by uncommenting the next line. O/W it is not worth
  % dealing with limited access restrictions and get/set
  % functions.
  % we would also need to deal with a no parameter constructor -- the
  % parameters would be something line zeros(3,1)
  % NM=class(NM,'NoiseModel');

