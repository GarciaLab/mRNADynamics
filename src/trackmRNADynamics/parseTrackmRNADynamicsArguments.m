function [app, retrack, optionalResults, displayFigures] =...
    parseTrackmRNADynamicsArguments(varargin)

app = {};
retrack = 1;
optionalResults = '';
displayFigures = 0;

for i = 2:length(varargin)
    if ~ischar(varargin{i})
        continue
    else
        if strcmpi(varargin{i}, 'app')
            app{1} = varargin{i + 1};
            app{2} = varargin{i + 2};
        elseif contains(varargin{i}, 'noRetrack', 'IgnoreCase', true)
            retrack = 0;
        elseif strcmpi(varargin{i}, 'optionalResults')
            optionalResults = varargin{i+1};
        elseif strcmpi(varargin{i}, 'displayFigures')
            displayFigures = 1;
        end
    end
    
    
end

end
