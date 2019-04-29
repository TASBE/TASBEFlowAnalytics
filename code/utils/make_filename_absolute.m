function absolute_name = make_filename_absolute(stem, path)

if is_octave()
    absolute = is_absolute_filename(stem);
else
    javaFileObj = javaObject('java.io.File', end_with_slash(stem));
    absolute = javaFileObj.isAbsolute();
end

if absolute
    absolute_name = end_with_slash(stem);
else
    absolute_name = end_with_slash(fullfile(path, stem));
end
