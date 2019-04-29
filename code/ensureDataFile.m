% ensureDataFile is a function in the DataFile class which returns a DataFile
% version of the input. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function DF = ensureDataFile(datafile)
   if isempty(datafile)
       DF = [];
       return 
   elseif isa(datafile,'DataFile')
       DF = datafile;
   elseif isa(datafile,'char')
       DF = DataFile(datafile);
   else
       TASBESession.error('TASBE:DataFile','InvalidType','Only files of char type can create DataFile objs');
   end
end
