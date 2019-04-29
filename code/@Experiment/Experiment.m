% EXPERIMENT holds experiment files and description. Can have only a
% subset of channels, but three parameters should be passed to the
% constructor (even if some are '') the constructor checks for
% parameters of the Channel class. Contains the following properties:
%         ExperimentName
%         InducerNames
%         InducerLevelsToFiles % This is an array with the first columns being inducer levels and the last column a list of filenames
%         % Below is not yet handled:
%         ProteinName % e.g. Tal
%         Construct % text for defining the construct that produced the data
%         Conditions % text for explaining the experiment conditions
%         Notes % text for anything related to this experiment
%         ExperimentDate
%         ExperiemntTime
%
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

 function E = Experiment(Name, InducerNames, InducerLevelsToFiles)
            E.ExperimentName=[];
            E.InducerNames=[];
            E.InducerLevelsToFiles=[]; 
            E.ProteinName =[]; 
            E.Construct =[];
        
            if nargin > 0
                E.ExperimentName = Name;
                E.InducerNames = InducerNames;
                E.InducerLevelsToFiles = InducerLevelsToFiles;
            end;
            
            E=class(E,'Experiment');
            
            % check to make sure InducerLevelsToFiles has the correct dimensions
            if size(E.InducerLevelsToFiles, 2) ~= 2 && size(E.InducerLevelsToFiles, 2) ~= 0
                TASBESession.error('TASBE:Experiment', 'DimensionMismatch', 'Transfer Curve analysis invoked with incorrect number of columns. Make sure InducerLevelsToFiles is a n X 2 matrix.');
            end  
           
            
% C1 = Channel('red','TxRed','r')
% example: E = Experiment({'Dox'}, {0, {'Dox0.fcs'}; 1, {'Dox1.cs'}; 10, {'Dox10.fcs', 'Dox10b.fcs'}}, C1, '', '')
