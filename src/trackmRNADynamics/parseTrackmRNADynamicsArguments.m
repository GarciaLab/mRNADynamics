function [app, retrack, optionalResults, displayFigures, trackByProximity] =...
    parseTrackmRNADynamicsArguments(varargin)

app = {};
retrack = false;
optionalResults = '';
displayFigures = false;
trackByProximity = false;

for i = 2:length(varargin)
    if ~ischar(varargin{i})
        continue
    else
        if strcmpi(varargin{i}, 'app')
            app{1} = varargin{i + 1};
            app{2} = varargin{i + 2};
        elseif strcmpi(varargin{i}, 'retrack')
            retrack = true;
        elseif strcmpi(varargin{i}, 'optionalResults')
            optionalResults = varargin{i+1};
        elseif strcmpi(varargin{i}, 'displayFigures')
            displayFigures = true;
        elseif strcmpi(varargin{i}, 'trackByProximity')
            trackByProximity = true;  
        end
    end
    
    
end

end
