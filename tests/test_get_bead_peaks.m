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
incorrect_batch = 'AP02';
vague_batch = 'Lot AD';
no_batch = '';

% test correct inputs
expected_peaks1 = [NaN, 692, 2192, 6027, 17491, 35673, 126914, 291016];
obtained_peaks1 = get_bead_peaks(correct_model, correct_channel, correct_batch1);
expected_peaks2 = [72, 646, 1705, 4828, 15995, 47602, 135873, 273002];
obtained_peaks2 = round(get_bead_peaks(correct_model, correct_channel, correct_batch2));
obtained_peaks3 = get_bead_peaks(correct_model, correct_channel, no_batch);

assertElementsAlmostEqual(expected_peaks1, obtained_peaks1);
assertElementsAlmostEqual(expected_peaks2, obtained_peaks2);
assertElementsAlmostEqual(expected_peaks1, obtained_peaks3);

% test incorrect model inputs
assertError(@()get_bead_peaks(incorrect_model, correct_channel, correct_batch1), 'get_bead_peaks:NoModel', 'No error was raised.');
assertError(@()get_bead_peaks(no_model, correct_channel, correct_batch1), 'get_bead_peaks:NoModel', 'No error was raised.');

% test incorrect channel inputs
assertError(@()get_bead_peaks(correct_model, incorrect_channel, correct_batch1), 'get_bead_peaks:NoChannel', 'No error was raised.');
assertError(@()get_bead_peaks(correct_model, no_channel, correct_batch1), 'get_bead_peaks:NoChannel', 'No error was raised.');

% test incorrect batch inputs
assertError(@()get_bead_peaks(correct_model, correct_channel, incorrect_batch), 'get_bead_peaks:NoBatch', 'No error was raised.');
assertError(@()get_bead_peaks(correct_model, correct_channel, vague_batch), 'get_bead_peaks:VagueInput', 'No error was raised.');
