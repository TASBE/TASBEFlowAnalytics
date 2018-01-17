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
            if isempty(log), log = {}; end;
            
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
    end
    
    methods(Static)
        function out = error(classname,name,message,varargin)
            event.name = name;
            event.classname = classname;
            event.type = 'error';
            event.message = sprintf(message,varargin{:});
            out = TASBESession.access('insert',classname,event);
            error([classname ':' name],event.message);
        end
        
        function out = warn(classname,name,message,varargin)
            event.name = name;
            event.classname = classname;
            event.type = 'failure';
            event.message = sprintf(message,varargin{:});
            out = TASBESession.access('insert',classname,event);
            warning([classname ':' name],event.message);
        end
        
        function out = skip(classname,name,message,varargin)
            event.name = name;
            event.classname = classname;
            event.type = 'skip';
            event.message = sprintf(message,varargin{:});
            out = TASBESession.access('insert',classname,event);
        end
        
        function out = succeed(classname,name,message,varargin)
            event.name = name;
            event.classname = classname;
            event.type = 'success';
            event.message = sprintf(message,varargin{:});
            out = TASBESession.access('insert',classname,event);
        end
        
        function reset()
            TASBESession.access('reset');
        end
        
        function out = list()
            out = TASBESession.access();
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
                fprintf(fid,out);
                fclose(fid);
            end
        end
    end
end
