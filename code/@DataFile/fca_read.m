function [fcsdat, fcshdr, fcsdatscaled] = fca_read(datafile)
% [fcsdat, fcshdr, fcsdatscaled] = fca_read(filename);
%
% Read CSV or FCS of flow cytometry data file and put the list mode  
% parameters to the fcsdat array with size of [NumOfPar TotalEvents]. 
% Some important header data are stored in the fcshdr structure:
% TotalEvents, NumOfPar, name

% If noarg was supplied
if nargin == 0
    TASBESession.warn('Read', 'NoDataFile', 'No DataFile obj provided! Returning empty dataset.');
    fcsdat = []; fcshdr = []; fcsdatscaled = [];
    return;
else
    if isempty(datafile.file)
        TASBESession.warn('Read','NoFile','No file provided! Returning empty dataset.'); 
        fcsdat = []; fcshdr = []; fcsdatscaled = [];
        return;
    end
    filecheck = dir(datafile.file);
    if size(filecheck,1) == 0
        TASBESession.warn('Read','NoFile',[datafile.file,': The file does not exist! Returning empty dataset.']); 
        fcsdat = []; fcshdr = []; fcsdatscaled = [];
        return;
    end
end

% If filename arg. only contain PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool

if strcmp(datafile.type, 'csv')
    [fcsdat, fcshdr, fcsdatscaled] = fca_readcsv(datafile.file, datafile.header, TASBEConfig.get('flow.maxEvents'));
else
    [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(datafile.file, TASBEConfig.get('flow.maxEvents'));
end

