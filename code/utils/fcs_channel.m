function c = fcs_channel()
% function c = fcs_channel()
%
% Makes a blank FCS file channel structure, of the sort produced by 
% fca_readfcs and fca_readcsv
%  
% Copyright (C) 2019-2019, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%  
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

c = struct();
c.name = '';
c.rawname = '';
c.range = [];
c.bit = [];
c.voltage = [];
c.filter = [];
c.gain = [];
c.emitter_wavelength = [];
c.emitter_power = [];
c.percent_light = [];
c.detector = [];
c.decade = [];
c.log = [];
c.logzero = [];
c.calibration = []; % FCS 3.1 field not yet properly parsed by fca_readfcs. Might unify with unit later
% Custom additions for harmonization with TASBE channels
c.unit = '';
c.print_name = '';
