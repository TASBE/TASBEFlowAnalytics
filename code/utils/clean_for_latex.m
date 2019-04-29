% CLEAN_FOR_LATEX performs cleanup on latex names (e.g., escaping underscores)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function cleaned = clean_for_latex(name)

cleaned = regexprep(name,'_','\\_');
