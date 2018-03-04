% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function outputfig(h,name,path)
% outputfig is an indirection that allows us to swap around the preferred
% output format with little effort

% Default path is current path
if nargin < 3, path = './'; end; % Note: frontslash works for both Windows and Mac/Unix
if path(numel(path)) ~= '/', path(numel(path)+1) = '/'; end; % ensure path ends in slash

if TASBEConfig.get('testing.fakeFigureSaves'), 
    warning('TASBE:Utilities','Test mode: not actually saving figure %s/%s',path,name);
    return;
end;

% If directory doesn't exist, try to create it
if ~isdir(path),
    warning('TASBE:Utilities','Directory does not exist, attempting to create it: %s',path);
    mkdir(path);
end

name = sanitize_filename(name);
%   Note: path must be added after sanitization, because path
%   characters need to be removed in sanitization

maxlength = 75; % 80 character database limit minus file extension
if (length(name) > maxlength)
    % shorten the name
    s = length(name)-maxlength+2;
    shortname = ['X', name(s:end)];
else
    shortname = name;
end

if (is_octave)
    % Prettify the fonts
    FN = findall(h,'-property','FontName');
    set(FN,'FontName','Helvetica');
    %FS = findall(h,'-property','FontSize');
    %set(FS,'FontSize',8);
    
    % print PDFs
    print(h,'-dpdfwrite',[path shortname]);
else
    print(h,'-depsc2',[path shortname]); % Still produce a vector format when invoked by hand via Matlab
end
print(h,'-dpng',[path shortname]);
saveas(h,[path shortname '.fig']);

if(strcmp(get(h,'visible'),'off'))
    close(h);
end
