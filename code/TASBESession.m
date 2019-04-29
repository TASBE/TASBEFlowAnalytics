% Copyright (C) 2011 - 2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE Flow Analytics distribution's top directory.
%
% This file is part of the TASBE Flow Analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the BBN Flow Cytometry
% package distribution's top directory.

classdef TASBESession
    methods(Static,Hidden)
        function out = access(mode,suite,entry)
            persistent log;
            if isempty(log), 
                log = {}; 
                % log initialization event
                suiteentry.name = 'TASBE:Session-1';
                suiteentry.contents = {};
                log{end+1} = suiteentry;
                event.name = 'Initialize';
                event.classname = 'TASBE:Session';
                try
                    version = tasbe_version();
                    event.message = sprintf('TASBE session logging enabled, release %i.%i.%s',version{1},version{2},version{3});
                    event.type = 'success';
                    log{end}.contents{end+1} = event;
                catch e
                    event.type = 'error';
                    event.message = 'Could not determine release version';
                    log{end}.contents{end+1} = event;
                    throw e;
                end
            end;
            
            % If there is no arguments, just return the whole thing for inspection
            if nargin==0, out = log; return; end
            % if the key is 'reset', then clear the log
            if strcmp(mode,'reset'), log = {}; return;
            elseif strcmp(mode,'insert'),
                suitename = sprintf('%s-%i',suite,numel(log));
                if(isempty(log) || ~strcmp(log{end}.name,suitename))
                    suiteentry.name = sprintf('%s-%i',suite,numel(log)+1);
                    suiteentry.contents = {};
                    log{end+1} = suiteentry;
                end
                log{end}.contents{end+1} = entry;
                out = entry;
            else
                error('TASBE:Session','Bad logging mode: %s',mode);
            end
        end
        
        function out = test_to_xml(event)
            if(strcmp(event.type,'success'))
                contents = sprintf('   <system-out>%s</system-out>\n',event.message);
            else
                contents = sprintf('   <%s>%s</%s>\n',event.type,event.message,event.type);
            end
            format = '  <testcase classname="%s" name="%s" time="0">\n%s  </testcase>\n';
            out = sprintf(format,event.classname,event.name,contents);
        end
        
        function out = test_to_xml_excel(event)
            if(strcmp(event.type,'success'))
                contents = sprintf('   <system-out>%s</system-out>\n',event.message);
            else
                contents = sprintf('   <%s>%s</%s>\n',event.type,event.message,event.type);
            end
            format = '  <testcase classname="%s" name="%s" time="%s">\n%s  </testcase>\n';
            out = sprintf(format,event.classname,event.name,datestr(now,'dd-mm-yyyy HH:MM:SS'),contents);
        end
        
        function out = suite_to_xml(suite)
            teststr = cell(numel(suite.contents),1);
            errs = 0; fails = 0;
            for i=1:numel(suite.contents)
                if strcmp(suite.contents{i}.type,'failure'), fails = fails+1;
                elseif strcmp(suite.contents{i}.type,'error'), errs = errs+1;
                end
                teststr{i} = TASBESession.test_to_xml(suite.contents{i});
            end
            tests = sprintf('%s',teststr{:});
            % TODO: should also include timestamp, compute time
            attributes = sprintf('errors="%i" tests="%i" failures="%i" time="0"',errs,numel(suite.contents),fails);
            out = sprintf(' <testsuite name="%s" %s>\n%s </testsuite>\n',suite.name,attributes,tests);
        end
        
        function out = suite_to_xml_excel(suite)
            teststr = cell(numel(suite.contents),1);
            errs = 0; fails = 0;
            for i=1:numel(suite.contents)
                if strcmp(suite.contents{i}.type,'failure'), fails = fails+1;
                elseif strcmp(suite.contents{i}.type,'error'), errs = errs+1;
                end
                teststr{i} = TASBESession.test_to_xml_excel(suite.contents{i});
            end
            tests = sprintf('%s',teststr{:});
            % TODO: should also include timestamp, compute time
            out = sprintf(' <testsuite>\n%s </testsuite>\n',tests);
        end

        function off = checkIfWarningOff(msgId)
            off = false;
            warnStruct = warning();
            for i=1:numel(warnStruct)
                if strcmp(msgId,warnStruct(i).identifier) && strcmp('off',warnStruct(i).state),
                    off = true;
                    return;
                end
            end
        end
    end
    
    methods(Static)
        function out = error(classname,name,message,varargin)
            % check if previous event is an abort, if so keep aborting
            log = TASBESession.list();
            if strcmp(log{end}.contents{end}.name, 'Abort')
                abort();
            else
                event.name = name;
                event.classname = classname;
                event.type = 'error';
                event.message = sprintf(message,varargin{:});
                out = TASBESession.access('insert',classname,event);
                errorStruct.message = strrep(event.message,'%','%%');
                errorStruct.identifier = [classname ':' name];
                error(errorStruct);
            end
        end
        
        function out = abort()
            event.name = 'Abort';
            event.classname = 'TASBESession';
            event.type = 'error';
            event.message = 'Exiting out of current analysis';
            out = TASBESession.access('insert',event.classname,event);
            errorStruct.message = strrep(event.message,'%','%%');
            errorStruct.identifier = [event.classname ':' event.name];
            error(errorStruct);
        end
        
        function out = warn(classname,name,message,varargin)
            % check if previous event is an abort, if so keep aborting
            log = TASBESession.list();
            if strcmp(log{end}.contents{end}.name, 'Abort')
                abort();
            else
                % abort if warning is turned off
                if TASBESession.checkIfWarningOff([classname ':' name]), return; end;
                % otherwise, continue
                event.name = name;
                event.classname = classname;
                event.type = 'failure';
                event.message = sprintf(message,varargin{:});
                out = TASBESession.access('insert',classname,event);
                warning([classname ':' name],strrep(event.message,'%','%%'));
            end
        end
        
        function out = skip(classname,name,message,varargin)
            % check if previous event is an abort, if so keep aborting
            log = TASBESession.list();
            if strcmp(log{end}.contents{end}.name, 'Abort')
                abort();
            else
                % abort if warning is turned off
                if TASBESession.checkIfWarningOff([classname ':' name]), return; end;
                % otherwise, continue
                event.name = name;
                event.classname = classname;
                event.type = 'skip';
                event.message = sprintf(message,varargin{:});
                out = TASBESession.access('insert',classname,event);
            end
        end
        
        function out = succeed(classname,name,message,varargin)
            % check if previous event is an abort, if so keep aborting
            log = TASBESession.list();
            if strcmp(log{end}.contents{end}.name, 'Abort')
                abort();
            else
                % abort if warning is turned off
                if TASBESession.checkIfWarningOff([classname ':' name]), return; end;
                % otherwise, continue
                event.name = name;
                event.classname = classname;
                event.type = 'success';
                event.message = sprintf(message,varargin{:});
                out = TASBESession.access('insert',classname,event);
                fprintf([strrep(event.message,'%','%%') '\n']);
            end
        end
        
        function out = notify(classname,name,message,varargin)
            % check if previous event is an abort, if so keep aborting
            log = TASBESession.list();
            if strcmp(log{end}.contents{end}.name, 'Abort')
                abort();
            else
                % abort if warning is turned off
                if TASBESession.checkIfWarningOff([classname ':' name]), return; end;
                % otherwise, continue
                event.name = name;
                event.classname = classname;
                event.type = 'success';
                event.message = sprintf(message,varargin{:});
                out = TASBESession.access('insert',classname,event);
                fprintf(['Note: ' strrep(event.message,'%','%%') '\n']);
            end
        end
        
        function reset()
            TASBESession.access('reset');
        end
        
        function out = list()
            out = TASBESession.access();
        end
        
        function selected = getLast(classname,name)
            warnings = fliplr(TASBESession.list());
            selected = [];
            for i=1:numel(warnings)
                suite = fliplr(warnings{i}.contents);
                for j=1:numel(suite)
                    event = suite{j};
                    if(strcmp(classname,event.classname) && strcmp(name,event.name))
                        selected = event;
                        return;
                    end
                end
            end
        end
        
        function out = to_xml(filename)
            contents = TASBESession.list();
            suitestr = cell(numel(contents),1);
            for i=1:numel(contents)
                suitestr{i} = TASBESession.suite_to_xml(contents{i});
            end
            suites = sprintf('%s',suitestr{:});
            header = '<?xml version="1.0" encoding="UTF-8"?>';
            out = sprintf('%s\n<testsuites>\n%s</testsuites>\n',header,suites);
            
            if nargin>0
                fid = fopen(filename,'w');
                fprintf(fid,strrep(out,'%','%%'));
                fclose(fid);
            end
        end
        
        function out = to_xml_excel(filename)
            contents = TASBESession.list();
            suitestr = cell(numel(contents),1);
            for i=1:numel(contents)
                suitestr{i} = TASBESession.suite_to_xml_excel(contents{i});
            end
            suites = sprintf('%s',suitestr{:});
            header = '<?xml version="1.0" encoding="UTF-8"?>';
            out = sprintf('%s\n<testsuites>\n%s</testsuites>\n',header,suites);
            
            if nargin>0
                fid = fopen(filename,'w');
                fprintf(fid,out);
                fclose(fid);
            end
        end
    end
end
