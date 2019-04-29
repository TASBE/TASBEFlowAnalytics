function test_suite = test_TASBEConfig
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_config

assert(TASBEConfig.isSet('foo') == false);

missingError = false;
try
    TASBEConfig.getexact('foo.bar.baz.qux')
    missingError = true;
catch e  % error is expected
end
if missingError, error('Should have failed on missing preference'); end;

assert(TASBEConfig.getexact('foo.bar.baz.qux',1) == 1);
log = TASBESession.list();
assertEqual(log{end}.contents{end}.name, 'UnknownSetting');
assert(TASBEConfig.getexact('foo.bar.baz.qux') == 1);

assert(TASBEConfig.getexact('foo.bar.baz.qux',2) == 1);
assert(TASBEConfig.get('foo.bar.baz.qux') == 1);

assert(TASBEConfig.set('foo.bar.baz.qux',3) == 3);
assert(TASBEConfig.get('foo.bar.baz.qux') == 3);

assert(TASBEConfig.list.foo.bar.baz.qux == 3);

try
    TASBEConfig.getexact('calibration.plotPath');
    missingError = true;
catch e  % error is expected
end
if missingError, error('Should have failed on missing preference'); end;

TASBEConfig.set('plots.plotPath','plots/');
assert(strcmp(TASBEConfig.get('calibration.plotPath'),'plots/'));

assert(TASBEConfig.isSet('foo') == true);
TASBEConfig.clear('foo');
assert(TASBEConfig.isSet('foo') == false);

try
    TASBEConfig.get('foo')
    missingError = true;
catch e  % error is expected
end
if missingError, error('Should have failed on missing preference'); end;

try
    TASBEConfig.list.foo.bar.baz.qux
    missingError = true;
catch e  % error is expected
end
if missingError, error('Should have failed on missing preference'); end;

%%%%
% Now test checkpointing
assert(TASBEConfig.isSet('foo') == false);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'test') && numel(l)==2 && strcmp(l{1},'test') && strcmp(l{2},'init'));

TASBEConfig.checkpoint('one');
TASBEConfig.set('foo',2);
assert(TASBEConfig.get('foo')==2);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'one') && numel(l)==3 && strcmp(l{1},'one') && strcmp(l{2},'test') && strcmp(l{3},'init'));


TASBEConfig.checkpoint('two');
TASBEConfig.set('foo',3);
assert(TASBEConfig.get('foo')==3);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'two') && numel(l)==4 && strcmp(l{1},'two') && strcmp(l{2},'one') && strcmp(l{3},'test') && strcmp(l{4},'init'));

TASBEConfig.checkpoint('two');
TASBEConfig.clear('foo');
assert(TASBEConfig.isSet('foo')==false);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'two') && numel(l)==4 && strcmp(l{1},'two') && strcmp(l{2},'one') && strcmp(l{3},'test') && strcmp(l{4},'init'));

TASBEConfig.checkpoint('two');
assert(TASBEConfig.get('foo')==2);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'two') && numel(l)==4 && strcmp(l{1},'two') && strcmp(l{2},'one') && strcmp(l{3},'test') && strcmp(l{4},'init'));

TASBEConfig.checkpoint('one');
assert(TASBEConfig.isSet('foo')==false);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'one') && numel(l)==3 && strcmp(l{1},'one') && strcmp(l{2},'test') && strcmp(l{3},'init'));

TASBEConfig.set('foo',5);
assert(TASBEConfig.get('foo')==5);

TASBEConfig.checkpoint('one');
assert(TASBEConfig.isSet('foo')==false);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'one') && numel(l)==3 && strcmp(l{1},'one') && strcmp(l{2},'test') && strcmp(l{3},'init'));

TASBEConfig.checkpoint('test');
assert(TASBEConfig.isSet('foo')==false);
[c, l] = TASBEConfig.checkpoints();
assert(strcmp(c,'test') && numel(l)==2 && strcmp(l{1},'test') && strcmp(l{2},'init'));


function test_config_serialization

output = TASBEConfig.to_json();
splitout = strsplit(output,'\n');

targets = {
    '"heatmapPlotType": "image"';
    '"beadChannel": "FITC",';
    };

for i=1:numel(targets),
    found = false;
    for j=1:numel(splitout)
        if(strcmp(strtrim(splitout{j}),targets{i})), 
            %fprintf('Found target: %s\n',targets{i});
            found = true; 
            break; 
        end;
    end
    assert(found);
end



function test_config_read

TASBEConfig.checkpoint('test');

% set up some scratch values
TASBEConfig.set('foo.bar',1);
TASBEConfig.set('foo.baz',2);
TASBEConfig.set('foo.qux',3);
TASBEConfig.set('bar.qux','quux');
TASBEConfig.set('baz','qux');
TASBEConfig.set('another',7);

% save the old values to JSON
old = TASBEConfig.to_json();

% load JSON to replace some of the scratch values
json = '{ "foo": { "bar": 4, "qux": 5 }, "bar": {"qux": 6}, "baz":"replaced"}';
TASBEConfig.load_from_json(json);
% confirm replacements
assertEqual(TASBEConfig.get('foo.bar'),4);
assertEqual(TASBEConfig.get('foo.baz'),2);
assertEqual(TASBEConfig.get('foo.qux'),5);
assertEqual(TASBEConfig.get('bar.qux'),6);
assertEqual(TASBEConfig.get('baz'),'replaced');
assertEqual(TASBEConfig.get('another'),7);

% reload old values from JSON and confirm:
TASBEConfig.load_from_json(old);
assertEqual(TASBEConfig.get('foo.bar'),1);
assertEqual(TASBEConfig.get('foo.baz'),2);
assertEqual(TASBEConfig.get('foo.qux'),3);
assertEqual(TASBEConfig.get('bar.qux'),'quux');
assertEqual(TASBEConfig.get('baz'),'qux');
assertEqual(TASBEConfig.get('another'),7);
