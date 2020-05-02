function performPostNuclearTracking(Prefix,...
    expandedAnaphaseFrames, nWorkers, schnitzcellsFile,...
    ellipsesFile, postTrackingSettings)
%
%This function performs additional operations on schnitzcells and Ellipses
%including 1) Calculating nuclear fluorescence, 2) Correcting tracking
%mistakes, 3) Adding useful relative times and positions for downstream
%analysis

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment);
Ellipses = getEllipses(liveExperiment);

if postTrackingSettings.fish
    schnitzcells = rmfield(schnitzcells, {'P', 'E', 'D'});
end

% Stitch the schnitzcells using Simon's code
if ~postTrackingSettings.noStitch
    disp('stitching schnitzes')
    [schnitzcells, Ellipses] = StitchSchnitz(Prefix, nWorkers);
end

%Extract the nuclear fluorescence values if we're in the right experiment
%type
if postTrackingSettings.intFlag
    schnitzcells = integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo);
end

for s = 1:length(schnitzcells)
    midFrame = ceil(length(schnitzcells(s).frames)/2);
    dif = double(schnitzcells(s).frames(midFrame)) - expandedAnaphaseFrames;
    cycle = find(dif>0, 1, 'last' );
    schnitzcells(s).cycle = uint8(cycle);
end

schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells,...
    FrameInfo, expandedAnaphaseFrames);

%perform some quality control
schnitzcells = filterSchnitz(schnitzcells,...
    [liveExperiment.yDim, liveExperiment.xDim]);

if ~exist([liveExperiment.resultsFolder,filesep,'APDetection.mat'], 'file')
    postTrackingSettings.shouldConvertToAP = false;
end

if postTrackingSettings.shouldConvertToAP
    [EllipsePos, APAngle, APLength]...
        = convertToFractionalEmbryoLength(Prefix);
    
    for s = 1:length(schnitzcells)
        for f = 1:length(schnitzcells(s).frames)
            ellipseInd = schnitzcells(s).cellno(f);
            schnitzcells(s).APPos(f) = EllipsePos{f}(ellipseInd);
        end
    end
end


save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


end