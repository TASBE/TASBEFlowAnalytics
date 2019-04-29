function test_suite = test_errors_excel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_errors_excel_filenames
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/faulty_batch_template1.xlsx', [end_with_slash(filepath) '../']);
    CM = load_or_make_testing_colormodel();
    assertExceptionThrown(@()batch_analysis_excel(extractor, CM), 'getExcelFilename:FilenameNotFound', 'No error was raised');
    
function test_errors_excel_conditions
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/faulty_batch_template1.xlsx', [end_with_slash(filepath) '../']);
    log = TASBESession.list();
    assertEqual(log{end}.contents{end}.name, 'InvalidValue');

function test_errors_excel_session
    if is_octave()
        assertExceptionThrown(@()analyzeFromExcel('', ''), 'Octave:invalid-index', 'No error was raised');   
    else
        assertExceptionThrown(@()analyzeFromExcel('', ''), 'MATLAB:badsubscript', 'No error was raised'); 
    end
    try
        assertExceptionThrown(@()analyzeFromExcel('test_templates/faulty_batch_template1.xlsx', 'comparative'), 'TASBE:Analysis:ColumnDimensionMismatch', 'No error was raised');
    catch
        assertExceptionThrown(@()analyzeFromExcel('test_templates/faulty_batch_template1.xlsx', 'comparative'), 'TASBE:ColorModel:MissingChannel', 'No error was raised');
    end
    assertExceptionThrown(@()analyzeFromExcel('test_templates/faulty_batch_template1.xlsx', 'none'), 'analyzeFromExcel:InvalidType', 'No error was raised');
    
function test_errors_excel_channel_roles
    CM = load_or_make_testing_colormodel();
    % Test 0 channel roles
    extractor = TemplateExtraction('test_templates/test_batch_template3.xlsx');
    assertExceptionThrown(@()getChannelRoles(CM, extractor), 'getChannelRoles:NoInformation', 'No error was raised'); 
    % Test 1 channel role
    extractor = TemplateExtraction('test_templates/test_batch_template2.xlsx');
    assertExceptionThrown(@()getChannelRoles(CM, extractor), 'getChannelRoles:TooFew', 'No error was raised'); 
    % Test 2 channel roles
    extractor = TemplateExtraction('test_templates/test_batch_template5.xlsx');
    [channel_roles, print_names] = getChannelRoles(CM, extractor); 
    expected_names = {'EYFP', 'mKate'};
    assertEqual(numel(channel_roles), 1);
    assertEqual(print_names, expected_names);
    % Test 3 channel roles
    extractor = TemplateExtraction('test_templates/test_batch_template1.xlsx');
    [channel_roles, print_names] = getChannelRoles(CM, extractor); 
    expected_names = {'EYFP', 'mKate', 'EBFP2'};
    assertEqual(numel(channel_roles), 1);
    assertEqual(print_names, expected_names);
    % Test 4 channel roles
    extractor = TemplateExtraction('test_templates/faulty_batch_template1.xlsx');
    [channel_roles, print_names] = getChannelRoles(CM, extractor); 
    expected_names = {'EYFP', 'mKate', 'EBFP2', 'EYFP'};
    assertEqual(numel(channel_roles), 2);
    assertEqual(print_names, expected_names); 
    
        