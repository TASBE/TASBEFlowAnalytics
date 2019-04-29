function directories_added=tasbe_set_path()
% sets the search path for TASBEFlowAnalytics
%
% directories_added=tasbe_set_path()
%
% Output:
%    directories_added          Cell string with directories that were
%                               added to the search path. If no
%                               directories were added to the search path
%                               (i.e., all functions for TASBEFlowAnalytics were
%                               already in the search path), then this is
%                               an empty cell array.
% Note:
%   - This adds directory in which this function resides, as well as the
%     'util' sub-directory to the search path

    return_directories_added=nargout>=1;

    root_dir=fileparts(mfilename('fullpath'));

    sub_dirs={'code', 'code/utils', 'library_code', 'library_code/jsonlab'}; % Add directories here as needed
    
    % Check if it's Octave; if not, assume matlab
    v = ver();
    isOctave = strcmpi(v(1).Name,'octave');
    if isOctave, 
        sub_dirs{end+1} = 'matlab_compat'; 
    else
        sub_dirs{end+1} = 'octave_compat'; 
    end;
    % Check if Matlab gmdistribution function is available; if not, load the free alternative
    if isempty(which('gmdistribution')), sub_dirs{end+1} = 'gmdistribution'; end;
    
    
    % Construct the path to be added by adding the pathsep() after
    % each subdirectory
    make_full_path_with_sep=@(sub_dir)[fullfile(root_dir,sub_dir) ...
                                        pathsep()];

    full_dirs=cellfun(make_full_path_with_sep,sub_dirs,...
                        'UniformOutput',false);

    if return_directories_added
        % store original path elements
        orig_path_cell=get_path_cell();
    end

    addpath(cat(2,full_dirs{:}));

    if return_directories_added
        new_path_cell=get_path_cell();
        directories_added=setdiff(new_path_cell, orig_path_cell);
    end

function pc=get_path_cell()
    pc=regexp(path(), pathsep(), 'split');

