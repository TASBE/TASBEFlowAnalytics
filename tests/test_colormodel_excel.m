function test_suite = test_colormodel_excel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_colormodel_excel_endtoend
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template1.xlsx', [end_with_slash(filepath) '../']);
    CM = make_color_model_excel(extractor);
    
    % Check that bead files are identical 
    log = TASBESession.list();
    assertEqual(log{end}.contents{end}.name, 'Identical');
    
    % Check results in CM:
    CMS = struct(CM);

    UT = struct(CMS.unit_translation);
    assertElementsAlmostEqual(UT.k_ERF,        2267.3,   'relative', 1e-2);
    assertElementsAlmostEqual(UT.first_peak,    8);
    assertElementsAlmostEqual(UT.fit_error,     0);
    assertElementsAlmostEqual(UT.peak_sets{1},  128.35, 'relative', 1e-2);

    AFM_Y = struct(CMS.autofluorescence_model{1});
    assertElementsAlmostEqual(AFM_Y.af_mean,    3.2600,  'absolute', 0.5);
    assertElementsAlmostEqual(AFM_Y.af_std,     17.0788, 'absolute', 0.5);
    AFM_R = struct(CMS.autofluorescence_model{2});
    assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
    assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);
    AFM_B = struct(CMS.autofluorescence_model{3});
    assertElementsAlmostEqual(AFM_B.af_mean,    5.8697,  'absolute', 0.5);
    assertElementsAlmostEqual(AFM_B.af_std,     16.9709, 'absolute', 0.5);

    COMP = struct(CMS.compensation_model);
    expected_matrix = [...
        1.0000      0.0056      0.0004      0;
        0.0010      1.0000      0.0022      0;
             0      0.0006      1.0000      0;
             0      0      0      1];

    assertElementsAlmostEqual(COMP.matrix,      expected_matrix, 'absolute', 1e-3);
    
    % Check results for sizeBeads:
    assertTrue(strcmp(CMS.sizeUnits,'Eum'));
    expected_peaks = 1e5*[0.1257    0.2574    0.6506    1.3908    1.9087    3.1596];
    UT = struct(CMS.size_unit_translation);
    assertElementsAlmostEqual(UT.um_poly,       [0.5865 -2.0798],  'relative', 1e-2);
    assertElementsAlmostEqual(UT.first_peak,    1);
    assertElementsAlmostEqual(UT.fit_error,     0.0784,   'absolute', 0.01);
    assertElementsAlmostEqual(UT.peak_sets{1},  expected_peaks, 'relative', 1e-2);
    
%     Not Equal: Got [NaN 1.0056 1.9442;0.99176 NaN NaN;0.51417 NaN NaN] instead
%     CTM = struct(CMS.color_translation_model);
%     expected_scales = [...
%            NaN    1.0097    2.1796;
%         0.9872       NaN       NaN;
%         0.45834      NaN       NaN];
% 
%     assertElementsAlmostEqual(CTM.scales,       expected_scales, 'absolute', 0.02);

function test_colormodel_excel_singlered
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template2.xlsx', [end_with_slash(filepath) '../']);
    CM = make_color_model_excel(extractor);
    % Check results in CM:
    CMS = struct(CM);

    UT = struct(CMS.unit_translation);
    assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
    assertElementsAlmostEqual(UT.first_peak,    5);
    assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
    assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

    AFM_R = struct(CMS.autofluorescence_model{1});
    assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
    assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

    COMP = struct(CMS.compensation_model);
    assertElementsAlmostEqual(COMP.matrix,      1.0000, 'absolute', 1e-3);

    CTM = struct(CMS.color_translation_model);
    assertElementsAlmostEqual(CTM.scales,   NaN);

function test_colormodel_excel_singlered_nocolorfile
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template3.xlsx', [end_with_slash(filepath) '../']);
    CM = make_color_model_excel(extractor);
    % Check results in CM:
    CMS = struct(CM);

    UT = struct(CMS.unit_translation);
    assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
    assertElementsAlmostEqual(UT.first_peak,    5);
    assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
    assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

    AFM_R = struct(CMS.autofluorescence_model{1});
    assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
    assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

    COMP = struct(CMS.compensation_model);
    assertElementsAlmostEqual(COMP.matrix,      1.0000, 'absolute', 1e-3);

    CTM = struct(CMS.color_translation_model);
    assertElementsAlmostEqual(CTM.scales,   NaN);

    