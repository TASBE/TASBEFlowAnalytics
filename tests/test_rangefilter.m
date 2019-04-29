function test_suite = test_rangefilter
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_range_filter

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';
blankfile = DataFile('fcs',[stem0312 'blank_P3.fcs']);
[~, hdr] = fca_read(blankfile);

hdr.par = hdr.par([7 10 11]); % 'FITC-A', 'PE-Tx-Red-YG-A', 'Pacific Blue-A'
data = [[2:10 1]; 2:2:20; 3:3:30]'; % some fake data to play with

% remove things outside of the permitted ranges
RF = RangeFilter('FITC-A',[2 9]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (2:9)'));
assert(all(Dnew(:,2) == (2:2:16)'));
assert(all(Dnew(:,3) == (3:3:24)'));

% should be a no-op from applying twice:
Dnewer = applyFilter(RF,hdr,Dnew);
assert(all(Dnewer(:,1) == (2:9)'));
assert(all(Dnewer(:,2) == (2:2:16)'));
assert(all(Dnewer(:,3) == (3:3:24)'));

% tighter range:
RF = RangeFilter('FITC-A',[4 5]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (4:5)'));
assert(all(Dnew(:,2) == (6:2:8)'));
assert(all(Dnew(:,3) == (9:3:12)'));

% open ranges:
RF = RangeFilter('PE-Tx-Red-YG-A',[-inf 5]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (2:3)'));
assert(all(Dnew(:,2) == (2:2:4)'));
assert(all(Dnew(:,3) == (3:3:6)'));

RF = RangeFilter('Pacific Blue-A',[10 inf]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == [(5:10) 1]'));
assert(all(Dnew(:,2) == (8:2:20)'));
assert(all(Dnew(:,3) == (12:3:30)'));

% two-component 'And'
RF = RangeFilter('Pacific Blue-A',[10 inf],'PE-Tx-Red-YG-A',[-inf 12]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (5:7)'));
assert(all(Dnew(:,2) == (8:2:12)'));
assert(all(Dnew(:,3) == (12:3:18)'));

% explicit inclusion of 'And'
RF = RangeFilter('Mode','And','Pacific Blue-A',[10 inf],'PE-Tx-Red-YG-A',[-inf 12]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (5:7)'));
assert(all(Dnew(:,2) == (8:2:12)'));
assert(all(Dnew(:,3) == (12:3:18)'));

RF = RangeFilter('Pacific Blue-A',[10 inf],'PE-Tx-Red-YG-A',[-inf 12],'Mode','And');
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == (5:7)'));
assert(all(Dnew(:,2) == (8:2:12)'));
assert(all(Dnew(:,3) == (12:3:18)'));

% use of 'Or' instead:
RF = RangeFilter('Pacific Blue-A',[22 inf],'Mode','Or','FITC-A',[-inf 3]);
Dnew = applyFilter(RF,hdr,data);
assert(all(Dnew(:,1) == [2 3 9 10 1]'));
assert(all(Dnew(:,2) == [2 4 16 18 20]'));
assert(all(Dnew(:,3) == [3 6 24 27 30]'));

% test application as part of a color model:
RF = RangeFilter('FITC-A',[10 inf]);
CM = load_or_make_testing_colormodel();
rawdat = read_filtered_au(CM,blankfile);
CM = add_prefilter(CM,RF);
filtdat = read_filtered_au(CM,blankfile);

assert(min(rawdat(:,7))<0);
assert(min(filtdat(:,7))>=10);
assert(min(rawdat(:,10))<0);
assert(min(filtdat(:,10))<0);

RF = RangeFilter('FITC-A',[1e5 inf]);
CM = load_or_make_testing_colormodel();
rawdat = readfcs_compensated_ERF(CM,blankfile,false,true);
CM = add_postfilter(CM,RF);
filtdat = readfcs_compensated_ERF(CM,blankfile,false,true);

assert(min(rawdat(:,1))<1e4);
assert(min(filtdat(:,1))>=1e5);
assert(min(rawdat(:,3))<1e4);
assert(min(filtdat(:,3))<1e4);
