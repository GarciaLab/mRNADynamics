function nNuclei = countNuclei(Prefix, varargin)

cycle = [];
displayFigures = false;

if length(varargin) == 1
    cycle = varargin{1};
elseif length(varargin) == 2
    displayFigures = true;
end

%just returns the number of nuclei in a cycle. 

%AR- i use this to verify that my labeling of the cycles is correct (e.g.
%~16 nuclei in my fov is definitely cycle 12)


[~,ProcPath,DropboxFolder,~, PreProcPath,~, Prefix, ~,~,~,~,~, ~, ~, movieDatabase]...
    = readMovieDatabase(Prefix);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase, DVResolution] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

ncFrames = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
diffFrames = [diff(ncFrames), NaN]; %cycle 14 not supported
midFrames = round(diffFrames/2) + ncFrames;


load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']);


nNuclei = [];
for nc = 9:13
    if midFrames(nc) > 0
        nNuclei(nc) = length(Ellipses{midFrames(nc)});
    end
end

nNuclei = [nNuclei, NaN]; %just to give it the dimensions 1-14

if ~isempty(cycle)
    nNuclei = nNuclei(cycle);
end

figure();
if displayFigures
    nNucleiAllFrames = [];
    for frame = 1:length(Ellipses)
        nNucleiAllFrames(frame) = length(Ellipses{frame});
    end
    plot(1:length(Ellipses), nNucleiAllFrames);
    for nc = 9:13
        if midFrames(nc) ~= 0
            xline(midFrames(nc), '-', ['Cycle ',num2str(nc)]);
        end   
    end
    xlabel('frame');
    ylabel('nuclei');
end



end