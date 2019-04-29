function CM = load_or_make_testing_colormodel2()

if exist('CM120312_2.mat','file')
    load('CM120312_2.mat');
    return;
end

%%%%%%%%%%%%%%%%
% if it doesn't already exist, make it:
TASBEConfig.set('plots.plotPath',[tempdir 'plots']);

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';


beadfile = [stem0312 'Beads_P3.fcs'];
blankfile = [stem0312 'blank_P3.fcs'];

% Autodetect gating with an N-dimensional gaussian-mixture-model
autogate = GMMGating(blankfile);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
% Channel takes FCS channel name, laser frequency (nm), filter center (nm), filter width (nm)
% Do not duplicate laser/filter information, as this may cause analysis collisions
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
% channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed
colorfiles{1} = [stem0312 'EYFP_P3.fcs']; % If there is only one channel, the color file is optional

channels{2} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{2} = setPrintName(channels{2}, 'mKate');
channels{2} = setLineSpec(channels{2}, 'r');
colorfiles{2} = [stem0312 'mkate_P3.fcs'];

channels{3} = Channel('Pacific Blue-A', 405, 450, 50);
channels{3} = setPrintName(channels{3}, 'EBFP2');
channels{3} = setLineSpec(channels{3}, 'b');
colorfiles{3} = [stem0312 'ebfp2_P3.fcs'];

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};
% Entries are: channel1, channel2, constitutive channel, filename
% This allows channel1 and channel2 to be converted into one another.
% If you only have two colors, you can set consitutive-channel to equal channel1 or channel2
colorpairfiles{1} = {channels{1}, channels{2}, channels{3}, [stem0312 'mkate_EBFP2_EYFP_P3.fcs']};
colorpairfiles{2} = {channels{1}, channels{3}, channels{2}, [stem0312 'mkate_EBFP2_EYFP_P3.fcs']};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
% Can also set bead channel if, for some reason, you don't want to use fluorescein as standard
% This defaults to FITC as it is strongly recommended to use fluorescein standards.
% TASBEConfig.set('beadChannel','FITC');

% Ignore all bead data below 10^[rangeMin] as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
%TASBEConfig.set('beads.peakThreshold', 200);
CM=set_ERF_channel_name(CM, 'FITC-A');
% Ignore channel data for ith channel if below 10^[value(i)]
TASBEConfig.set('colortranslation.channelMinimum',[2,2,2]);

% When dealing with very strong fluorescence, use secondary channel to segment
% TASBEConfig.set('beads.secondaryBeadChannel','PE-Tx-Red-YG-A');
CM = add_prefilter(CM,autogate);

% Execute and save the model
CM=resolve(CM);
save('-V7','CM120312_2.mat','CM');
