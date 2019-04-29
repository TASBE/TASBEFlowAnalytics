function [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(filename,clip_events)
% [fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(filename);
%
% Read FCS 2.0, 3.0 and 3.1 type flow cytometry data file and put the list mode  
% parameters to the fcsdat array with the size of [NumOfPar TotalEvents]. 
% Some important header data are stored in the fcshdr structure:
% TotalEvents, NumOfPar, starttime, stoptime and specific info for parameters
% as name, range, bitdepth, logscale(yes-no) and number of decades.
%
% [fcsdat, fcshdr] = fca_readfcs;
% Without filename input the user can select the desired file(s)
% using the standard open file dialog box.
%
% [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename);
% Supplying the third output the fcsdatscaled array contains the scaled     
% parameters. It might be useful for logscaled parameters, but no effect 
% in the case of linear parameters. The log scaling is the following
% operation for the "ith" parameter:  
% fcsdatscaled(:,i) = ...
%   10.^(fcsdat(:,i)/fcshdr.par(i).range*fcshdr.par(i).decade;);
%
%[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(filename);
% In that case the script calculates the compensated fluorescence 
% intensities (fcsdat_comp) as well, if spillover data exist in the header 
%
% Version: 15/Nov/2017
% University of Debrecen, Institute of Medical Imaging
% Laszlo Balkay 
% balkay.laszlo@med.unideb.hu
%
% History
% 14/08/2006 I made some changes in the code by the suggestion of 
% Brian Harms <brianharms@hotmail.com> and Ivan Cao-Berg <icaoberg@cmu.edu> 
% (given at the user reviews area of Mathwork File exchage) The program should work 
% in the case of Becton EPics DLM FCS2.0, CyAn Summit FCS3.0 and FACSDiva type 
% list mode files.
%
% 29/01/2008 Updated to read the BD LSR II file format and including the comments of
% Allan Moser (Cira Discovery Sciences, Inc.)
%
% 24/01/2009 Updated to read the Partec CyFlow format file. Thanks for
% Gavin A Price (GAP)
% 
% 20/09/2010 Updated to read the Accuri C6 format file. Thanks for
% Rob Egbert, University of Washington
%
% 07/11/2011 Updated to read Luminex 100 data file. Thanks for
% Ofir Goldberger, Stanford University
%
% 11/05/2013 The fluorescence compensation is implemeted into the code.
% Thanks for Rick Stanton, J. Craig Venter Institute, La Jolla, San Diego
%
% 12/02/2013  
%   MOre accurate compensation correction and amplification gain scaling is added. 
% Thanks for Rachel Finck(RLF); Garry Nolan's lab at Stanford University
%   Appropriate byte offset for the data segment is included for large 
%   file size (100Mbyte>).  
% Thanks for Andrea Pagnani(AP) /Politecnico Torino, Human Genetics Foundation
% and RLF
%
% 16/05/2014
% The FCS 3.0 standard enables the mixture of word lengths in the data, this
% upgrade modified the code according to. The linefeed (ASCII code 10) as 
% the mnemonic separator was also added. 
% Thanks for William Peria /Fred Hutchinson Cancer Research Center
%
% 15/04/2017
% First try for the FCS 3.1 standard. The files were supported by 
% Matteo Ferla, Department of Biochemistry, University of Oxford.
%
% 1/11/2017
% Minor update to allow FCS3.1 data from Attune NxT to be imported
% Not sure if it is fully FCS3.1 compliant, 
% Thanks for Robin Cleveland, Biomedical Engineering, Oxford
% 
% Further updated by Jacob Beal, 2010 - 2018



% if noarg was supplied
if nargin == 0
     [FileName, FilePath] = uigetfile('*.*','Select fcs2.0 file');
     filename = [FilePath,FileName];
     if FileName == 0;
          fcsdat = []; fcshdr = []; fcsdatscaled = []; fcsdat_comp= [];
          return;
     end
else
    if isempty(filename)
        TASBESession.warn('FCS:Read','NoFile','No file provided! Returning empty dataset.'); 
        fcsdat = []; fcshdr = []; fcsdatscaled = [];
        return;
    end
    % Remove the NULL ascii character. This makes strange things!
    filename_asciicode = int16(filename);
    filename(find(filename_asciicode==0)) = [];
    % check if the file exists
    filecheck = dir(filename);
    if size(filecheck,1) == 0 || size(filecheck,1) >1
        TASBESession.warn('FCS:Read','NoFile',[filename,': The file does not exist! Returning empty dataset.']); 
        fcsdat = []; fcshdr = []; fcsdatscaled = []; fcsdat_comp= [];
        return;
    end
end
if nargin<2, clip_events = 1e6; end;

% if filename arg. only contain PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool
[FilePath, FileNameMain, fext] = fileparts(filename);
FilePath = [FilePath filesep];
FileName = [FileNameMain, fext];
if  isempty(FileNameMain)
    currend_dir = cd;
    cd(FilePath);
    [FileName, FilePath] = uigetfile('*.*','Select FCS file');
     filename = [FilePath,FileName];
     if FileName == 0;
          fcsdat = []; fcshdr = []; fcsdatscaled= []; fcsdat_comp= [];
          return;
     end
     cd(currend_dir);
end

%fid = fopen(filename,'r','ieee-be');
fid = fopen(filename,'r','b');
fcsheader_1stline   = fread(fid,64,'char');
fcsheader_type = char(fcsheader_1stline(1:6)');

%% reading the header
if strcmp(fcsheader_type,'FCS1.0')
    TASBESession.error('FCS:Read','FCS1',[FileName ': FCS 1.0 file type is not supported.']);
    fcsdat = []; fcshdr = []; fcsdatscaled= []; fcsdat_comp= [];
    fclose(fid);
    return;
elseif  strcmp(fcsheader_type,'FCS2.0') || strcmp(fcsheader_type,'FCS3.0') || strcmp(fcsheader_type,'FCS3.1') % FCS2.0 or FCS3.0 or FCS3.1 types
    fcshdr.fcstype = fcsheader_type;
    FcsHeaderStartPos   = str2num(char(fcsheader_1stline(11:18)'));
    FcsHeaderStopPos    = str2num(char(fcsheader_1stline(19:26)'));
    FcsDataStartPos     = str2num(char(fcsheader_1stline(27:34)'));

    status = fseek(fid,0,'bof');
    fcsheader_total = fread(fid,FcsHeaderStopPos+1,'char');%read the total header
    status = fseek(fid,FcsHeaderStartPos,'bof');
    fcsheader_main = fread(fid,FcsHeaderStopPos-FcsHeaderStartPos+1,'char');%read the main header
    warning off MATLAB:nonIntegerTruncatedInConversionToChar;
    fcshdr.filename = FileName;
    fcshdr.filepath = FilePath;
    % "The first character of the primary TEXT segment contains the
    % delimiter" (FCS standard)
    if fcsheader_main(1) == 12
        mnemonic_separator = 'FF';
    elseif fcsheader_main(1) == 9 
        mnemonic_separator = 'TAB'; 
    elseif fcsheader_main(1) == 10
        mnemonic_separator = 'LF';
    elseif fcsheader_main(1) == 30
        % "Record Selector" asci code = 30
        rs_pos = find(fcsheader_total == 30);
        fcsheader_total(rs_pos)  =  int8('|');
        rs_pos = find(fcsheader_main == 30);
        fcsheader_main(rs_pos)  =  int8('|');
        mnemonic_separator = '|';
    else
        mnemonic_separator = char(fcsheader_main(1));
    end
    %
    % if the file size larger than ~100Mbyte the previously defined
    % FcsDataStartPos = 0. In that case the $BEGINDATA parameter stores the correct value
    %
    if ~FcsDataStartPos
         FcsDataStartPos = str2num(get_mnemonic_value('$BEGINDATA',fcsheader_main, mnemonic_separator));
    end
    %
    if mnemonic_separator == '@';% WinMDI
        TASBESession.error('FCS:Read','WinMDI',[FileName ': The file can not be read (Unsupported FCS type: WinMDI histogram file)']);
        fcsdat = []; fcshdr = [];fcsdatscaled= []; fcsdat_comp= [];
        fclose(fid);
        return;
    end
    fcshdr.TotalEvents = str2num(get_mnemonic_value('$TOT',fcsheader_main, mnemonic_separator));
    if fcshdr.TotalEvents == 0
        fcsdat = 0;
        fcsdatscaled = 0;
        return
    end
    fcshdr.NumOfPar = str2num(get_mnemonic_value('$PAR',fcsheader_main, mnemonic_separator));
%     if strcmp(mnemonic_separator,'LF')
%         fcshdr.NumOfPar = fcshdr.NumOfPar + 1;
%     end
%     
    fcshdr.Creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    
    % comp matrix reader
    comp = get_mnemonic_value('SPILL',fcsheader_main,mnemonic_separator); 
    if ~isempty(comp)
        compcell=regexp(comp,',','split');
        nc=str2double(compcell{1});        
        fcshdr.CompLabels=compcell(2:nc+1);
        fcshdr.CompMat=reshape(str2double(compcell(nc+2:end)'),[nc nc])';       
    else
        fcshdr.CompLabels=[];
        fcshdr.CompMat=[]; 
    end
    plate = get_mnemonic_value('PLATE NAME',fcsheader_main,mnemonic_separator);
    if ~isempty(plate)
        fcshdr.plate=plate;
    end
    
    for i=1:fcshdr.NumOfPar
        fcshdr.par(i) = fcs_channel();
        fcshdr.par(i).unit = 'a.u.'; % FCS file units are always a.u. until calibrated
        fcshdr.par(i).name = get_mnemonic_value(['$P',num2str(i),'N'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).rawname = fcshdr.par(i).name;
        fcshdr.par(i).range = str2num(get_mnemonic_value(['$P',num2str(i),'R'],fcsheader_main, mnemonic_separator));
        fcshdr.par(i).bit = str2num(get_mnemonic_value(['$P',num2str(i),'B'],fcsheader_main, mnemonic_separator));

        %%%%
        % Pick up additional (optional) parameters: - JSB
        %%%%
        % $PnV Detector voltage for parameter n.
        fcshdr.par(i).voltage = get_mnemonic_value(['$P',num2str(i),'V'],fcsheader_main, mnemonic_separator);
        if(fcshdr.par(i).voltage), fcshdr.par(i).voltage = str2num(fcshdr.par(i).voltage); end;
        % $PnF Name of optical filter for parameter n.
        fcshdr.par(i).filter = get_mnemonic_value(['$P',num2str(i),'F'],fcsheader_main, mnemonic_separator);
        % $PnG Amplifier gain used for acquisitionof parameter n.
        fcshdr.par(i).gain = get_mnemonic_value(['$P',num2str(i),'G'],fcsheader_main, mnemonic_separator);
        if(fcshdr.par(i).gain), fcshdr.par(i).gain = str2num(fcshdr.par(i).gain); end;
        % $PnL Excitation wavelength for parameter n.
        fcshdr.par(i).emitter_wavelength = get_mnemonic_value(['$P',num2str(i),'L'],fcsheader_main, mnemonic_separator);
        if(fcshdr.par(i).emitter_wavelength), fcshdr.par(i).emitter_wavelength = str2num(fcshdr.par(i).emitter_wavelength); end;
        % $PnO Excitation power for parameter n.
        fcshdr.par(i).emitter_power = get_mnemonic_value(['$P',num2str(i),'O'],fcsheader_main, mnemonic_separator);
        if(fcshdr.par(i).emitter_power), fcshdr.par(i).emitter_power = str2num(fcshdr.par(i).emitter_power); end;
        % $PnP Percent of emitted light collectedby parameter n.
        fcshdr.par(i).percent_light = get_mnemonic_value(['$P',num2str(i),'P'],fcsheader_main, mnemonic_separator);
        if(fcshdr.par(i).percent_light), fcshdr.par(i).percent_light = str2num(fcshdr.par(i).percent_light); end;
        % $PnT Detector type for parameter n.
        fcshdr.par(i).detector = get_mnemonic_value(['$P',num2str(i),'T'],fcsheader_main, mnemonic_separator);
        %%% not using these:
        % $Pkn Peak channel number of univariatehistogram for parameter n.
        % $PKNn Count in peak channel of univariatehistogram for parameter n.
        % $PnS Name used for parameter n.
        replacename = get_mnemonic_value(['$P',num2str(i),'S'],fcsheader_main, mnemonic_separator);
        % Replace name if nickname is present and not just whitespace
        if(numel(replacename)>0 && ~isempty(find(~isspace(replacename), 1))), fcshdr.par(i).name = replacename; end;
        % FCS 3.1-only parameters:
        if(strcmp(fcsheader_type,'FCS3.1')), 
            fcshdr.par(i).calibration = get_mnemonic_value(['$P',num2str(i),'CALIBRATION'],fcsheader_main, mnemonic_separator);
        end

        %==============   Changed way that amplification type is treated 
        par_exponent_str= (get_mnemonic_value(['$P',num2str(i),'E'],fcsheader_main, mnemonic_separator));
        if isempty(par_exponent_str)
            % There is no "$PiE" mnemonic in the Lysys format
            % in that case the PiDISPLAY mnem. shows the LOG or LIN definition
            islogpar = get_mnemonic_value(['P',num2str(i),'DISPLAY'],fcsheader_main, mnemonic_separator);
            if strcmp(islogpar,'LOG')
               par_exponent_str = '5,1'; 
            else % islogpar = LIN case
                par_exponent_str = '0,0';
            end
        end
       
        par_exponent= str2num(par_exponent_str);
        fcshdr.par(i).decade = par_exponent(1);
        if fcshdr.par(i).decade == 0
            fcshdr.par(i).log = 0;
            fcshdr.par(i).logzero = 0;
        else
            fcshdr.par(i).log = 1;
            if (par_exponent(2) == 0)
              fcshdr.par(i).logzero = 1;
            else
              fcshdr.par(i).logzero = par_exponent(2);
            end
        end
        gain_str = get_mnemonic_value(['$P',num2str(i),'G'],fcsheader_main, mnemonic_separator); % added by RLF 
        if ~isempty(gain_str) 
            fcshdr.par(i).gain=str2double(gain_str);
        else
            fcshdr.par(i).gain=1;
        end
  
    end
    
    fcshdr.starttime = get_mnemonic_value('$BTIM',fcsheader_main, mnemonic_separator);
    fcshdr.stoptime = get_mnemonic_value('$ETIM',fcsheader_main, mnemonic_separator);
    fcshdr.timestep = get_mnemonic_value('$TIMESTEP',fcsheader_main, mnemonic_separator);
    fcshdr.cytometry = get_mnemonic_value('$CYT',fcsheader_main, mnemonic_separator);
    fcshdr.date = get_mnemonic_value('$DATE',fcsheader_main, mnemonic_separator);
    fcshdr.byteorder = get_mnemonic_value('$BYTEORD',fcsheader_main, mnemonic_separator);
    if strcmp(fcshdr.byteorder, '1,2,3,4')
        machineformat = 'ieee-le';
    elseif strcmp(fcshdr.byteorder, '4,3,2,1')
        machineformat = 'ieee-be';
    end
    fcshdr.datatype = get_mnemonic_value('$DATATYPE',fcsheader_main, mnemonic_separator);
    fcshdr.system = get_mnemonic_value('$SYS',fcsheader_main, mnemonic_separator);
    fcshdr.project = get_mnemonic_value('$PROJ',fcsheader_main, mnemonic_separator);
    fcshdr.experiment = get_mnemonic_value('$EXP',fcsheader_main, mnemonic_separator);
    fcshdr.cells = get_mnemonic_value('$Cells',fcsheader_main, mnemonic_separator);
    fcshdr.creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
else
    TASBESession.error('FCS:Read','UnsupportedFCSType',[FileName,': The file can not be read (Unsupported FCS type)']);
    fcsdat = []; fcshdr = []; fcsdatscaled= []; fcsdat_comp= [];
    fclose(fid);
    return;
end

% deal with size overflow bug for very large FCS files, in which data start gets incorrectly set to zero
% This appears to be triggered, with FACSdiva at least, when the FCS file is more than 100MB, which makes
% the data stop position more than 8 characters and causes it to collide with the start field.   -JSB
if FcsDataStartPos==0,
    TASBESession.warn('FCS:Read','BadDataStart',[FileName,': FCS file has invalid data start position; guessing based on header end']);
    FcsDataStartPos = FcsHeaderStopPos+6;
end

% optionally truncate events to avoid memory problems with extremely large FCS files -JSB
if fcshdr.TotalEvents>clip_events,
    TASBESession.warn('FCS:Read','TooManyEvents',[FileName,': FCS file has more than %i events (%i events); truncating to avoid memory problems'],clip_events,fcshdr.TotalEvents);
    fcshdr.TotalEvents = clip_events;
end

%% reading the events
status = fseek(fid,FcsDataStartPos,'bof');
if strcmp(fcsheader_type,'FCS2.0')
    if strcmp(mnemonic_separator,'\') || strcmp(mnemonic_separator,'FF')... %ordinary or FacsDIVA FCS2.0 
           || strcmp(mnemonic_separator,'/') || strcmp(mnemonic_separator,'TAB') % added by GAP 1/22/09 %added by RLF 09/02/10
        if fcshdr.par(1).bit == 16
            fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',machineformat)');
            fcsdat_orig = uint16(fcsdat);%//
            if strcmp(fcshdr.byteorder,'1,2')...% this is the Cytomics data
                    || strcmp(fcshdr.byteorder, '1,2,3,4') %added by GAP 1/22/09
                fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
            end
            
            new_xrange = 1024;
            for i=1:fcshdr.NumOfPar
                if fcshdr.par(i).range > 4096
                    fcsdat(:,i) = fcsdat(:,i)*new_xrange/fcshdr.par(i).range;
                    fcshdr.par(i).range = new_xrange;
                end
            end
        elseif fcshdr.par(1).bit == 32
                if fcshdr.datatype ~= 'F'
                    fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32')');
                else % 'LYSYS' case
                    fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32')');
                end
        else 
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end
    elseif strcmp(mnemonic_separator,'!');% Becton EPics DLM FCS2.0
        fcsdat_ = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16', 'ieee-le')';
        fcsdat = zeros(fcshdr.TotalEvents,fcshdr.NumOfPar);
        for i=1:fcshdr.NumOfPar
            bintmp = dec2bin(fcsdat_(:,i));
            fcsdat(:,i) = bin2dec(bintmp(:,7:16)); % only the first 10bit is valid for the parameter  
        end
    end
    fclose(fid);
elseif strcmp(fcsheader_type,'FCS3.0') || strcmp(fcsheader_type,'FCS3.1') 
    if strcmp(mnemonic_separator,'|') && strcmp(fcshdr.datatype,'I') && (strcmp( fcshdr.cytometry,'CyAn') || strcmp( fcshdr.cytometry,'MoFlo Astrios' ))% CyAn Summit FCS3.0
        fcsdat_ = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],['uint',num2str(fcshdr.par(1).bit)],machineformat)'; 
        %fcsdat_ = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',machineformat)');
        fcsdat = zeros(size(fcsdat_));
        new_xrange = 2^16; 
        for i=1:fcshdr.NumOfPar
            fcsdat(:,i) = fcsdat_(:,i)*new_xrange/fcshdr.par(i).range;
            fcshdr.par(i).range = new_xrange;
        end
    elseif strcmp(mnemonic_separator,'/')
        if findstr(lower(fcshdr.cytometry),'accuri')  % Accuri C6, this condition added by Rob Egbert, University of Washington 9/17/2010
            fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'int32',machineformat)');
        elseif findstr(lower(fcshdr.cytometry),'partec')%this block added by GAP 6/1/09 for Partec, copy/paste from above
            fcsdat = uint16(fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',machineformat)');
            %fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
        elseif findstr(lower(fcshdr.cytometry),'lx') % Luminex data
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'int32',machineformat)';
            fcsdat = mod(fcsdat,1024);
        else  %Assume FCS3.1 format that is the same as FCS3.0
            if strcmp(fcshdr.datatype,'D')
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'double',machineformat)';
            elseif strcmp(fcshdr.datatype,'F')
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',machineformat)';
            elseif strcmp(fcshdr.datatype,'I')
                fcsdat = fread(fid,[sum([fcshdr.par.bit]/16) fcshdr.TotalEvents],'uint16',machineformat)'; % sum: William Peria, 16/05/2014
            end
        end
    else % ordinary FCS 3.0
        if strcmp(fcshdr.datatype,'D')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'double',machineformat)';
        elseif strcmp(fcshdr.datatype,'F')
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',machineformat)';
        elseif strcmp(fcshdr.datatype,'I')
            if sum([fcshdr.par.bit])/fcshdr.par(1).bit ~= fcshdr.NumOfPar % if the bitdepth different at different pars
                fcsdat = fread(fid,[sum([fcshdr.par.bit]/16) fcshdr.TotalEvents],'uint16',machineformat)'; % sum: William Peria, 16/05/2014 
            else
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],['uint',num2str(fcshdr.par(1).bit)],machineformat)';
            end
        else 
            TASBESession.error('FCS:Read','UnsupportedDataType',[FileName,': Unsupported FCS 3.0 data type: ' fcshdr.datatype]);
        end
    end
    fclose(fid);
end

%% this is converting Partec to FacsDIVA_FCS20 format if save_FacsDIVA_FCS20 = 1;
save_FacsDIVA_FCS20 = 0;
if strcmp(fcshdr.cytometry ,'partec PAS') && save_FacsDIVA_FCS20    
    fcsheader_main2 = fcsheader_main;
    sep_place = strfind(char(fcsheader_main'),'/');
    fcsheader_main2(sep_place) = 12;
    fcsheader_1stline2 = fcsheader_1stline;
    fcsheader_1stline2(31:34) = double(num2str(FcsHeaderStopPos+1));
    fcsheader_1stline2(43:50) = double('       0');
    fcsheader_1stline2(51:58) = double('       0');
    FileSize =  length(fcsheader_main2(:))+ length(fcsheader_1stline2(1:FcsHeaderStartPos))+ 2*length(fcsdat_orig(:));
    space_char(1:8-length(num2str(FileSize)))= ' ';
    fcsheader_1stline2(35:42) = double([space_char,num2str(FileSize)]);

    fid2 = fopen([FilePath, FileNameMain,'_', fext],'w','b');
    fwrite(fid2,[fcsheader_1stline2(1:FcsHeaderStartPos)],'char');
    fwrite(fid2,fcsheader_main2,'char');
    fwrite(fid2,fcsdat_orig','uint16');
    fclose(fid2);
end

%% this is for Gy Vamosi converting MoFlo Astrios with changing some hdr parameters if save_MoFlo_Astrios = 1
save_MoFlo_Astrios = 0;
if strcmp(fcshdr.cytometry ,'MoFlo Astrios') && save_MoFlo_Astrios    
    fcsheader_total_char = char(fcsheader_total)';
    intmax_pos = strfind(fcsheader_total_char,'4294967296');
    for i=1:length(intmax_pos)
            fcsheader_total_char(intmax_pos(i):intmax_pos(i)+length('4294967296')-1) = ['65536|    '];
    end
%     byteorder_pos = strfind(fcsheader_total_char,'1,2,3,4');
%     fcsheader_total_char(byteorder_pos:byteorder_pos+length('1,2,3,4')-1) = '4,3,2,1';
    
%     datatype_pos = strfind(fcsheader_total_char,'DATATYPE|I');
%     fcsheader_total_char(datatype_pos:datatype_pos+length('DATATYPE|I')-1) = 'DATATYPE|F';
    
%     pexp_pos = strfind(fcsheader_total_char,'0.0,0.0');
%     for i=1:length(pexp_pos)
%         fcsheader_total_char(pexp_pos(i):pexp_pos(i)+length('0.0,0.0')-1) = '0,0|$BL';
%     end
    
%     pgain_pos = strfind(fcsheader_total_char,'G|1.0000000000');
%     for i=1:length(pgain_pos)
%         fcsheader_total_char(pgain_pos(i):pgain_pos(i)+length('G|1.0000000000')-1) = 'G|1.0|$BLLLLLL';
%     end
    
    fcsheader_total2 = int8(fcsheader_total_char');
    
    fid2 = fopen([FilePath, FileNameMain,'_', fext],'w','b');
    fwrite(fid2,fcsheader_total2,'char');
    fwrite(fid2,uint32(fcsdat'),'uint32','ieee-le');
    %fwrite(fid2,fcsdat','float32','ieee-le');
    fclose(fid2);
end

% Ensure FCS data is ordinary floating point numbers
fcsdat = double(fcsdat);

%%  calculate the scaled events (for log scales) 

fcsdatscaled = zeros(size(fcsdat));
if numel(fcsdat)>0 && nargout>2
    for  i = 1 : fcshdr.NumOfPar
        Xlogdecade = fcshdr.par(i).decade;
        XChannelMax = fcshdr.par(i).range;
        Xlogvalatzero = fcshdr.par(i).logzero;
        if fcshdr.par(i).gain~=1
            fcsdatscaled(:,i)  = double(fcsdat(:,i))./fcshdr.par(i).gain;
            
        elseif fcshdr.par(i).log
            fcsdatscaled(:,i) = Xlogvalatzero*10.^(double(fcsdat(:,i))/XChannelMax*Xlogdecade);
        else
            fcsdatscaled(:,i)  = fcsdat(:,i);
        end
    end
    
end

%% calculate the compensated events
if nargout>3 && ~isempty(fcshdr.CompLabels) 
    
    compcols=zeros(1,nc);
    colLabels={fcshdr.par.name};
    for i=1:nc
        compcols(i)=find(strcmp(fcshdr.CompLabels{i},colLabels));
    end
    fcsdatcomp = fcsdatscaled;
    fcsdatcomp(:,compcols) = fcsdatcomp(:,compcols)/fcshdr.CompMat;
else fcsdatcomp=[];
end



function mneval = get_mnemonic_value(mnemonic_name,fcsheader,mnemonic_separator)

if strcmp(mnemonic_separator,'\')  || strcmp(mnemonic_separator,'!') ...
        || strcmp(mnemonic_separator,'|') || strcmp(mnemonic_separator,'@')...
        || strcmp(mnemonic_separator, '/') ... % added by GAP 1/22/08
        || strcmp(mnemonic_separator, '*') % added by JSB 1/7/19. TODO: generalize to any delimiter
    mnemonic_startpos = findstr(char(fcsheader'),[mnemonic_name,mnemonic_separator]);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length;
    next_slashes = findstr(char(fcsheader(mnemonic_stoppos+1:end)'),mnemonic_separator);
    next_slash = next_slashes(1) + mnemonic_stoppos;

    mneval = char(fcsheader(mnemonic_stoppos+1:next_slash-1)');
elseif strcmp(mnemonic_separator,'FF')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 12);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;

    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
elseif strcmp(mnemonic_separator,'TAB')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 9);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
    
elseif strcmp(mnemonic_separator, 'LF')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_linefeeds = find( fcsheader(mnemonic_stoppos+1:end) == 10);
    next_linefeed = next_linefeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_linefeed-1)');
    
elseif strcmp(mnemonic_separator, 'RS')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_linefeeds = find( fcsheader(mnemonic_stoppos+1:end) == 30);
    next_linefeed = next_linefeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_linefeed-1)');

  
end
