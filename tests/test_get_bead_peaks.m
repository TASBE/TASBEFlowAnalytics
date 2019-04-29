function test_suite = test_get_bead_peaks
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
    
function test_beadPeaksEndtoend
correct_model = 'SpheroTech RCP-30-5A';
incorrect_model = 'Unknown';
no_model = '';

correct_channel = 'FITC';
incorrect_channel = 'TE';
no_channel = '';

correct_batch1 = 'Lot AA';
correct_batch2 = 'AG01';
correct_batch3 = 'Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R';
incorrect_batch = 'AP02';
vague_batch = 'Lot AD';
no_batch = '';

% test correct inputs
expected_peaks1 = [NaN, 692, 2192, 6027, 17491, 35673, 126914, 291016];
obtained_peaks1 = get_bead_peaks(correct_model, correct_channel, correct_batch1);
expected_peaks2 = [72, 646, 1705, 4828, 15995, 47602, 135873, 273002];
obtained_peaks2 = round(get_bead_peaks(correct_model, correct_channel, correct_batch2));
obtained_peaks3 = get_bead_peaks(correct_model, correct_channel, no_batch);
obtained_peaks4 = get_bead_peaks(correct_model, correct_channel, correct_batch3);

assertElementsAlmostEqual(expected_peaks1, obtained_peaks1);
assertElementsAlmostEqual(expected_peaks2, obtained_peaks2);
assertElementsAlmostEqual(expected_peaks1, obtained_peaks3);
assertElementsAlmostEqual(expected_peaks1, obtained_peaks4);

% test incorrect model inputs
assertExceptionThrown(@()get_bead_peaks(incorrect_model, correct_channel, correct_batch1), 'TASBE:BeadCatalog:NoModel', 'No error was raised.');
assertExceptionThrown(@()get_bead_peaks(no_model, correct_channel, correct_batch1), 'TASBE:BeadCatalog:NoModel', 'No error was raised.');

% test incorrect channel inputs
assertExceptionThrown(@()get_bead_peaks(correct_model, incorrect_channel, correct_batch1), 'TASBE:BeadCatalog:NoChannel', 'No error was raised.');
assertExceptionThrown(@()get_bead_peaks(correct_model, no_channel, correct_batch1), 'TASBE:BeadCatalog:NoChannel', 'No error was raised.');

% test incorrect batch inputs
assertExceptionThrown(@()get_bead_peaks(correct_model, correct_channel, incorrect_batch), 'TASBE:BeadCatalog:NoBatch', 'No error was raised.');
assertExceptionThrown(@()get_bead_peaks(correct_model, correct_channel, vague_batch), 'TASBE:BeadCatalog:VagueInput', 'No error was raised.');

%%%%%%%
% An all-combinations test with URCP-38-2K
correct_model2 = 'SpheroTech URCP-38-2K';
correct_batch4 = 'Lot AJ02';
correct_batch5 = 'Lot AJ03';
correct_batch6 = 'AJ02';
correct_batch7 = 'AJ03';

expected_peaks4567 = [190	3868	33701	97304	251992	485706];
obtained_peaks4 = get_bead_peaks(correct_model2, correct_channel, correct_batch4);
obtained_peaks5 = get_bead_peaks(correct_model2, correct_channel, correct_batch5);
obtained_peaks6 = get_bead_peaks(correct_model2, correct_channel, correct_batch6);
obtained_peaks7 = get_bead_peaks(correct_model2, correct_channel, correct_batch7);

assertElementsAlmostEqual(expected_peaks4567, obtained_peaks4,'relative',1e-2);
assertElementsAlmostEqual(expected_peaks4567, obtained_peaks5,'relative',1e-2);
assertElementsAlmostEqual(expected_peaks4567, obtained_peaks6,'relative',1e-2);
assertElementsAlmostEqual(expected_peaks4567, obtained_peaks7,'relative',1e-2);
