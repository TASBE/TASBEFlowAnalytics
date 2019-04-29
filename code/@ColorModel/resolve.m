% RESOLVE is called after the construction and configuration of a
% ColorModel object. It includes unit translation from beads and 
% autofluorescence, compensation, and color translation model generation.   
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM=resolve(CM)
    % fill in channel descriptors from designated file (default = beadfile)
    if TASBEConfig.isSet('flow.channel_template_file'), 
        template = ensureDataFile(TASBEConfig.get('flow.channel_template_file'));
    else template = CM.BeadFile; 
    end;
    [fcsdat fcshdr] = fca_read(template);
    % Remember channel descriptions, for later confirmation
    for i=1:numel(CM.Channels),
        [~, desc] = get_fcs_color(fcsdat,fcshdr,getName(CM.Channels{i}));
        CM.Channels{i} = setDescription(CM.Channels{i},desc);
        % TODO: figure out how to add FSC and SSC channel descriptions (used by filters) for confirmation
    end
    
    
    % build model
    % First, unit translation from beads
    if TASBEConfig.isSet('calibration.overrideUnits')
        k_ERF = TASBEConfig.get('calibration.overrideUnits');
        CM.unit_translation = UnitTranslation('Specified',k_ERF,[],[],{});
        CM.standardUnits = 'arbitrary units';
        TASBESession.warn('TASBE:ColorModel','OverrideBeads','Overriding units with specified k_ERF value of %d',k_ERF);
    else
        [UT, CM] = beads_to_ERF_model(CM,CM.BeadFile);
        CM.unit_translation = UT;
    end
    
    % Next, unit translation from size beads, if there are any
    if isempty(CM.um_channel_name)
        CM.size_unit_translation = [];
        CM.sizeUnits = 'no size channel defined';
        TASBESession.warn('TASBE:ColorModel','NoSizeChannel','No size bead calibration configured');
    else
        [UT, CM] = beads_to_um_model(CM,CM.SizeBeadFile);
        CM.size_unit_translation = UT;
    end
    
    % Next, autofluorescence and compensation model
    if TASBEConfig.isSet('calibration.overrideAutofluorescence')
        afmean = TASBEConfig.get('calibration.overrideAutofluorescence');
        if numel(afmean)==1, afmean = afmean*ones(numel(CM.Channels),1); end;
        TASBESession.warn('TASBE:ColorModel','OverrideAutofluorescence','Overriding autofluorescence model with specified values.');
        for i=1:numel(afmean),
            CM.autofluorescence_model{i} = AutoFluorescenceModel(CM.Channels{i},afmean(i)*ones(10,1));
            if(CM.Channels{i} == CM.ERF_channel)
                CM.autofluorescence_model{i}=ERFize(CM.autofluorescence_model{i},1,getK_ERF(CM.unit_translation));
            end
        end
    else
        CM.autofluorescence_model = computeAutoFluorescence(CM);
    end
    if TASBEConfig.isSet('calibration.overrideCompensation')
        matrix = TASBEConfig.get('calibration.overrideCompensation');
        TASBESession.warn('TASBE:ColorModel','OverrideCompensation','Overriding compensation model with specified values.');
        CM.compensation_model = LinearCompensationModel(matrix, zeros(size(matrix)));
    else
        CM.compensation_model = computeColorCompensation(CM);
    end
    CM.initialized = 0.5; % enough to read in AU
    if TASBEConfig.get('compensation.plot'), plot_compensated_controls(CM); end;
    
    % finally, color translation model
    if TASBEConfig.isSet('calibration.overrideTranslation')
        scales = TASBEConfig.get('calibration.overrideTranslation');
        color_translation_model = ColorTranslationModel(CM.Channels,scales);
        for i=1:numel(CM.Channels),
            if(CM.Channels{i}==CM.ERF_channel) i_ERF = i; end;
        end
        for i=1:numel(CM.Channels),
            if(CM.Channels{i}==CM.ERF_channel || isUnprocessed(CM.Channels{i})) continue; end;
            AFMi = CM.autofluorescence_model{i};
            k_ERF=getK_ERF(CM.unit_translation);
            CM.autofluorescence_model{i}=ERFize(AFMi,scales(i,i_ERF),k_ERF);
        end
        TASBESession.warn('TASBE:ColorModel','OverrideTranslation','Overriding translation scaling with specified values.');
    else
        [color_translation_model CM] = computeColorTranslations(CM);
    end
    CM.color_translation_model = color_translation_model;
    CM.initialized = 1; % enough to read in (pseudo)ERF
    
    [ok CM] = confirm_ERF_translations(CM);% warn if can't translate all to ERF
    
    CM.noise_model = computeNoiseModel(CM);
    
