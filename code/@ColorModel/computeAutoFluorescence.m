% COMPUTEAUTOFLUORESCENCE returns an array of AutoFluorescenceModel objects
% obtained by processing channels on the ColorModel's blank control.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function autofluorescence_model=computeAutoFluorescence(CM)

makePlots = TASBEConfig.get('autofluorescence.plot');

% Import data from FCS files
[rawfcs, fcshdr] = read_filtered_au(CM,CM.BlankFile);

autofluorescence_model = cell(numel(CM.Channels),1);
for i=1:numel(CM.Channels)
    channel = CM.Channels{i};
    name=getName(channel);
    if(isUnprocessed(channel))
        TASBESession.notify('TASBE:Autofluorescence','UnprocessedChannel','Skipping autofluorescence computation for unprocessed channel %s',name);
        continue;
    end
    found = false;
    for j=1:numel(fcshdr.par)
        if(strcmp(name,fcshdr.par(j).name) || strcmp(name,fcshdr.par(j).rawname))
            autofluorescence_model{i} = AutoFluorescenceModel(channel,rawfcs(:,j));
            if(channel == CM.ERF_channel) % ERFize model if possible
                autofluorescence_model{i}=ERFize(autofluorescence_model{i},1,getK_ERF(CM.unit_translation));
            end
            % optional plot
            if makePlots, plot_autofluorescence_control(autofluorescence_model{i},rawfcs(:,j)); end
            
            found = true; break;
        end
    end
    if(found==false)
        TASBESession.error('TASBE:Autofluorescence','MissingChannel','Unable to find required channel %s',name)
    end
end
