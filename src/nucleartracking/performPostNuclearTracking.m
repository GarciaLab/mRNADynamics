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

%the cellno field of schnitzcells is hugely problematic. let's delete it 
%and rebuild it with addschnitzindextoellipses
schnitzcells = rmfield(schnitzcells, 'cellno'); 


%we'll make sure cellnos and ellipses correspond well.
[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


% Stitch the schnitzcells using Simon's fantastic and clever code
if ~postTrackingSettings.noStitch
    disp('stitching schnitzes')
    [schnitzcells, Ellipses] = StitchSchnitzv2(Prefix, nWorkers);
end

%making copies for validation later on
ellipsesOld = Ellipses;
schnitzcellsOld = schnitzcells;


%add the length field to schnitzcells
for schnitzIndex = 1:length(schnitzcells)
    framesPerSchnitz = length(schnitzcells(schnitzIndex).frames);
    for frameIndex = 1:framesPerSchnitz
        radius = single(mean(Ellipses{schnitzcells(schnitzIndex).frames(frameIndex)}(...
            schnitzcells(schnitzIndex).cellno(frameIndex),3:4)));
        if ~isreal(radius)
            radius = nan;
            warning('non real radii returned for schnitz. not sure what happened here.');
        end
        schnitzcells(schnitzIndex).len = radius*ones(1, framesPerSchnitz);
    end
end


ellipsesSizeUnchanged(ellipsesOld, Ellipses); %just a test, not sure if it works
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

%Extract the nuclear fluorescence values if we're in the right experiment
%type
if postTrackingSettings.intFlag
    schnitzcells = integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo);
end


ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

for schnitzIndex = 1:length(schnitzcells)
    midFrame = ceil(length(schnitzcells(schnitzIndex).frames)/2);
    dif = double(schnitzcells(schnitzIndex).frames(midFrame)) - expandedAnaphaseFrames;
    cycle = find(dif>0, 1, 'last' );
    schnitzcells(schnitzIndex).cycle = uint8(cycle);
end



ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells,...
    FrameInfo, expandedAnaphaseFrames);



ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

%perform some quality control
schnitzcells = filterSchnitz(schnitzcells,...
    [liveExperiment.yDim, liveExperiment.xDim]);



ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

if ~exist([liveExperiment.resultsFolder,filesep,'APDetection.mat'], 'file')
    postTrackingSettings.shouldConvertToAP = false;
end

if postTrackingSettings.shouldConvertToAP
    [EllipsePos, APAngle, APLength]...
        = convertToFractionalEmbryoLength(Prefix);
    
    
    
    ellipsesSizeUnchanged(ellipsesOld, Ellipses);
    schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);
    
    for schnitzIndex = 1:length(schnitzcells)
        for frameIndex = 1:length(schnitzcells(schnitzIndex).frames)
            frame = schnitzcells(schnitzIndex).frames(frameIndex);
            ellipseInd = schnitzcells(schnitzIndex).cellno(frameIndex);
            schnitzcells(schnitzIndex).APPos(frameIndex) = EllipsePos{frame}(ellipseInd);
        end
    end
end


ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


end