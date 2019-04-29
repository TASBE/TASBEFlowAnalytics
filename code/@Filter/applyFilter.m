% APPLYFILTER is an abstract filter method filters raw fcs data using an
% inputted Filter object
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function gated = applyFilter(F,fcshdr,rawfcs)
% F is the filter object
% fcshdr is the header for the FCS file
% rawfcs is the data to be filtered
% gated is the replacement for rawfcs

TASBESEssion.error('TASBE:Filter','AbstractFilterInvoked','Attempted to invoke abstract filter method');
