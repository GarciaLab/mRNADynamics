function performPostNuclearTracking(Prefix,...
    expandedAnaphaseFrames, nWorkers, schnitzcellsFile,...
    ellipsesFile, postTrackingSettings)
%
%This function performs additional operations on schnitzcells and Ellipses
%including 1) Calculating nuclear fluorescence, 2) Correcting tracking
%mistakes, 3) Adding useful relative times and positions for downstream
%analysis

disp('Starting post-tracking procedures...')

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment);
Ellipses = getEllipses(liveExperiment);

if postTrackingSettings.fish
    schnitzcells = rmfield(schnitzcells, {'P', 'E', 'D'});
end

%%
%Optionally, tighten the nuclear contours and convert the circles to
%true ellipses
if postTrackingSettings.doAdjustNuclearContours
    if postTrackingSettings.chooseHis
        [hisFile, hisPath] = uigetfile([liveExperiment.procFolder, filesep,'*.*']);
        hisMat = imreadStack([hisPath, filesep, hisFile]);
    else
        hisMat = getHisMat(liveExperiment);
    end
    pixelSize_um = liveExperiment.pixelSize_um;
    
    disp('Adjusting nuclear contours...');

    parfor frame = 1:length(Ellipses)
%             for frame = 1:length(Ellipses)
        
        Ellipses{frame} = adjustNuclearContours(Ellipses{frame},...
            hisMat(:, :, frame), pixelSize_um);
        
    end
    
    disp('Nuclear contours adjusted.');

    %TrackNuclei handles empty frames poorly, so let's fill them in.
    Ellipses = fillEmptyXYFrames(Ellipses);
    
    save2(ellipsesFile, Ellipses);
    TrackNuclei(Prefix, 'nWorkers', nWorkers, 'retrack')
    return;
    
end
%%

%we'll make sure cellnos and ellipses correspond well.
% [Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);

[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);

[Ellipses, schnitzcells] = addStrayEllipsesToSchnitzcells(Ellipses, schnitzcells);


save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


% Stitch the schnitzcells using Simon's fantastic and clever code
if ~postTrackingSettings.noStitch
    [schnitzcells, Ellipses] = StitchSchnitzv3(Prefix, nWorkers);
end

[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);

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
if ~isempty(liveExperiment.inputChannels)
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

% %perform some quality control
% schnitzcells = filterSchnitz(schnitzcells,...
%     [liveExperiment.yDim, liveExperiment.xDim]);



ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);


try
    if exist([liveExperiment.resultsFolder,filesep,'APDetection.mat'], 'file')
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
catch
    warning('Failed to update nuclear AP positions for unknown reason.')
end


ellipsesSizeUnchanged(ellipsesOld, Ellipses);
schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);


end