function [sortByFrame, sortByLength, ForCompileAll, SpeedMode, SisterMode, ...
    ncRange, projectionMode, plot3DGauss, intScale, NC, ...
    startNC, endNC, optionalResults, nWorkers, fish] = determineCheckParticleTrackingOptions(varargin)
%DETERMINEOPTIONS Summary of this function goes here
%   Detailed explanation goes here


%Flag to sort or not particles according to their starting frame
sortByFrame=true;
%Flag to sort by the number of times a spot was found in each particle
sortByLength=false;
%Flag to just save the data. This is good for CompileAll
ForCompileAll=false;
%Flag to plot only ellipses for current particle & save time
SpeedMode = false;
%Decide whether you want to do sister chromatid analysis
SisterMode = false;
%Decide whether you want to only see nc13
ncRange = 0;
% This is for the projection mode
projectionMode = 'None';
%plot 3D gaussian fitting intensities in tracefig
plot3DGauss = false;
%optional results if you have multiple prefixes with different results
%folders
optionalResults = '';
nWorkers = 1;
fish = false;

intScale = 1;

% these variables are meaningless if ncRange is 0
NC = -1;
startNC = -1;
endNC = -1;

for i=1:length(varargin)
    if strcmpi(varargin{i},'NoSort')
        sortByFrame=0;
    elseif strcmpi(varargin{i},'sortByLength')
        sortByLength=1;
    elseif strcmpi(varargin{i},'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i},'ForCompileAll')
        ForCompileAll=1;
    elseif strcmpi(varargin{i}, 'speedmode')
        SpeedMode = 1;
    elseif strcmpi(varargin{i}, 'plot3DGauss')
        plot3DGauss = 1;
    elseif strcmpi(varargin{i}, 'sistermode')
        SisterMode = 1;
    elseif strcmpi(varargin{i}, 'intScale')
        intScale = varargin{i+1};
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
        elseif strcmpi(varargin{i}, 'fish')
        fish = true;
    elseif strcmpi(varargin{i},'nc') % checking for the desired nc range
        ncRange = 1;
        NC = varargin{i+1};
        % startNC and endNC will be the variable names that have the start of the nc(s) of interest
        if length((varargin{i+1})) == 2
            startNC = ['nc' num2str(varargin{i+1}(1))];
            endNC = ['nc' num2str(varargin{i+1}(2) +1)];% Not including the next nc
        else
            startNC = ['nc' num2str(varargin{i+1})];
            endNC = ['nc' num2str(varargin{i+1} + 1)]; % Not including the next nc
        end
    end
end


end

