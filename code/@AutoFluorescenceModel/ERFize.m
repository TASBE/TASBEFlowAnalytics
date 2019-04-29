% ERFize changes an AutoFluorescenceModel to be able to apply to 
% calibrated (ERF) values as well as arbitrary unit values
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
    
function AFM=ERFize(AFM,scale,k_ERF)
    AFM.af_mean_ERF = AFM.af_mean*scale*k_ERF;
    AFM.af_std_ERF = AFM.af_std*scale*k_ERF;
