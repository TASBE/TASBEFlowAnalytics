function test_suite = test_get_bead_models
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
    
function test_beadModelsEndtoend

expected_models = {'SpheroTech RCP-30-5A'
    'SpheroTech RCP-30-5'
    'Clontech AcGFP 632594'
    'Clontech mCherry 632595'
    'SpheroTech URCP-38-2K'};

results = get_bead_models();

assertEqual(expected_models{1}, results{1});
assertEqual(expected_models{2}, results{2});
assertEqual(expected_models{3}, results{3});
assertEqual(expected_models{4}, results{4});
assertEqual(expected_models{5}, results{5});