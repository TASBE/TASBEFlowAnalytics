% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function sanitized = sanitize_filename(name)
% remove all non-"word" characters from a name (e.g., a filename)

sanitized = regexprep(name,'[^\w+\-*]','');
if numel(sanitized) < numel(name)
    warning('TASBE:SanitizeFilename','Name "%s" contains unsafe characters: abbreviated to "%s"',name,sanitized);
end
