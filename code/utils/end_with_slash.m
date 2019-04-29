% ensure path ends in slash and slahes are forward not backward
function new_path = end_with_slash(path)
    path = strrep(path, '\', '/');
    if path(numel(path)) ~= '/'
        path(numel(path)+1) = '/'; 
    end 
    new_path = path;
end