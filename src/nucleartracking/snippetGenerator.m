function performPostNuclearTracking(Prefix, fish,...
    expandedAnaphaseFrames, nWorkers, shouldConvertToAP, schnitzcellsFile,...
    ellipsesFile)
%
%This function performs additional operations on schnitzcells and Ellipses
%including 1) Calculating nuclear fluorescence, 2) Correcting tracking
%mistakes, 3) Adding useful relative times and positions for downstream
%analysis

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment);
Ellipses = getEllipses(liveExperiment);

nFrames = length(Ellipses);


%Extract the nuclear fluorescence values if we're in the right experiment
%type
if intFlag
    schnitzcells = integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo);
end

if fish schnitzcells = rmfield(schnitzcells, {'P', 'E', 'D'}); end

if track && ~noBreak
    [schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(...
        schnitzcells, Ellipses, expandedAnaphaseFrames, nFrames);
    save2(ellipsesFile, Ellipses);
    save2(schnitzcellsFile, schnitzcells);
end

% Stitch the schnitzcells using Simon's code
if ~noStitch
    disp('stitching schnitzes')
    StitchSchnitz(Prefix, nWorkers);
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

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);

% try
Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
if shouldConvertToAP
    [EllipsePos, APAngle, APLength]...
        = convertToFractionalEmbryoLength(Prefix);
end
for s = 1:length(schnitzcells)
    for f = 1:length(schnitzcells(s).frames)
        ellipseInd = schnitzcells(s).cellno(f);
        schnitzcells(s).APPos(f) = EllipsePos{f}(ellipseInd);
    end
end
% end


save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


end