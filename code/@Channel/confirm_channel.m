% CONFIRM_CHANNEL checks to make sure an inputted channel matches with
% required description (returns boolean) in terms of parameters such as
% voltage, gain, and wavelength.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function ok = confirm_channel(channel,alt_desc)

ok = true;
desc = channel.description;

% numerical parameters:
if(numel(desc.voltage) && numel(alt_desc.voltage))
    if desc.voltage ~= alt_desc.voltage
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %d, but is %d',...
            desc.name,'voltage',desc.voltage,alt_desc.voltage);
        ok = false;
    end
end

if(numel(desc.gain) && numel(alt_desc.gain))
    if desc.gain ~= alt_desc.gain
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %d, but is %d',...
            desc.name,'gain',desc.gain,alt_desc.gain);
        ok = false;
    end
end

if(numel(desc.emitter_wavelength) && numel(alt_desc.emitter_wavelength))
    if desc.emitter_wavelength ~= alt_desc.emitter_wavelength
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %d, but is %d',...
            desc.name,'emitter wavelength',desc.emitter_wavelength,alt_desc.emitter_wavelength);
        ok = false;
    end
end

if(numel(desc.emitter_power) && numel(alt_desc.emitter_power))
    if desc.emitter_power ~= alt_desc.emitter_power
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %d, but is %d',...
            desc.name,'emitter power',desc.emitter_power,alt_desc.emitter_power);
        ok = false;
    end
end

if(numel(desc.percent_light) && numel(alt_desc.percent_light))
    if desc.percent_light ~= alt_desc.percent_light
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %d, but is %d',...
            desc.name,'percent light',desc.percent_light,alt_desc.percent_light);
        ok = false;
    end
end

% string parameters
if(numel(desc.filter) && numel(alt_desc.filter))
    if ~strcmp(desc.filter,alt_desc.filter)
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %s, but is %s',...
            desc.name,'filter',desc.filter,alt_desc.filter);
        ok = false;
    end
end

if(numel(desc.detector) && numel(alt_desc.detector))
    if ~strcmp(desc.detector,alt_desc.detector)
        TASBESession.warn('Channel:Description','ValueMismatch','In channel %s, descriptor ''%s'' expected to be %s, but is %s',...
            desc.name,'detector',desc.detector,alt_desc.detector);
        ok = false;
    end
end


