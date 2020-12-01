clear
close all

% get data structures
Prefix = '2020-03-13-2xDl-Ven_snaBAC-mCh_Leica_Zoom2_7uW14uW_09';
liveExperiment = LiveExperiment(Prefix);
Particles = getParticles(liveExperiment);
trackingOptions = getTrackingOptions(liveExperiment);
Spots = getSpots(liveExperiment);

% set basic parameters for auto-gap filling 
minTraceLength = 10;
nbSize = 7;
minDensity = 5/7;
maxDistance = 2;
nbFilter = ones(1,nbSize);

% track elligible frames for fitting, and corresponding position info
autoFitInfo = struct(...
                'Frame',[], ...                
                'xPos', [], ...
                'yPos', [], ...
                'zPos', []);

% iterate through particles and identify candidate frames/positions              
for Channel = 1:length(Particles)
    for p = 1:length(Particles{Channel})
        FrameVecFull = Particles{Channel}(p).framesFull;
        ObservedFrameFilter = Particles{Channel}(p).obsFrameFilter;
        
        if sum(ObservedFrameFilter) >= minTraceLength
            %%% Find points to be integrated
            % find islands of dense points
            islandPointVec = conv(ObservedFrameFilter,nbFilter,'same')/nbSize >= minDensity & ObservedFrameFilter;
            % find points that are a part of these islands
            memberPointVec = conv(islandPointVec,nbFilter,'same') >= 0 & ObservedFrameFilter;
            % find points that are elligible to be auto-integrated
            elligiblePointVec = bwdist(memberPointVec)<=maxDistance & ~ObservedFrameFilter;
            % get frames
            elligibleFrames = FrameVecFull(elligiblePointVec);
            % record info
            autoFitInfo(Channel).Frame = [autoFitInfo.Frame elligibleFrames];
            autoFitInfo(Channel).xPos = [autoFitInfo.xPos Particles{Channel}(p).xPosInf(elligiblePointVec)'];
            autoFitInfo(Channel).yPos = [autoFitInfo.yPos Particles{Channel}(p).yPosInf(elligiblePointVec)'];
            autoFitInfo(Channel).zPos = [autoFitInfo.zPos Particles{Channel}(p).zPosInf(elligiblePointVec)];            
        end
                
    end
end

% call auto-fitting function...