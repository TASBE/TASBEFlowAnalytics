% Constuctor for Channel class with properties:
% name
% Laser
% FilterSpec
% facs name, filter and laser which will be used for equivalance
% checking
% 'pseudo' marker is used for internal use when the controls were
% not correct, but the data is being used anyway (e.g., for early
% trials of a method
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function C = Channel(name, laser, filterC, filterW, pseudo,unprocessed)
    if nargin == 0 
        C.name ='';
        C.Laser = 0;
        C.FilterCenter = 0;
        C.FilterWidth = 0;
        C.PseudoUnits = 0;
        C.Units = 'a.u.';
    elseif nargin >= 4
        C.name = name;
        C.Laser = laser;
        C.FilterCenter = filterC;
        C.FilterWidth = filterW;
        if nargin >= 5
            C.PseudoUnits = pseudo;
        else
            C.PseudoUnits = 0;
        end
        C.Units = 'a.u.';
        % Warn if we're seeing unspecified channels
        if((nargin<6 || ~unprocessed) && (C.Laser==0 || C.FilterCenter==0 || C.FilterWidth==0)),
            TASBESession.warn('TASBE:Channel','Underspecified','Channel %s has unspecified laser and/or filter.  Unspecified channels may be confused together.',C.name);
        end
    end
    
    C.description = []; % file descriptor: will get filled in by colormodel resolution
    C.LineSpec='';
    C.PrintName='';
    if nargin<6, 
        C.unprocessed = false; 
        if(strcmp(C.name,'FSC') || strcmp(C.name,'FSC-A') || strcmp(C.name,'FSC-H') || strcmp(C.name,'FSC-W') || ...
            strcmp(C.name,'SSC') || strcmp(C.name,'SSC-A') || strcmp(C.name,'SSC-H') || strcmp(C.name,'SSC-W'))
            TASBESession.notify('TASBE:Channel','UnprocessedChannel','Channel %s has been automatically detected as unprocessed from its name',C.name);
            C.unprocessed=true;
        end
    else
        C.unprocessed=unprocessed; % set to true for FSC/SSC channels
    end
    
    C=class(C,'Channel');


