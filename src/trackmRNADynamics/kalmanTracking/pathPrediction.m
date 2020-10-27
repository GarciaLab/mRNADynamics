function ParticlesTemp = pathPrediction(ParticlesTemp, backwardTracks, trackingInfo, kalmanOptions)

  
    [frameVec, frameOrder] = sort(backwardTracks.Frame);
    
    % generate sorted position info
    framesFull = max([1,frameVec(1)-trackingInfo.nExtrapFrames]):min([trackingInfo.nFrames,frameVec(end)+trackingInfo.nExtrapFrames]);
    posData = NaN(length(framesFull), size(backwardTracks.MeasurementVec,2));
    posData(ismember(framesFull,frameVec),:) = backwardTracks.MeasurementVec(frameOrder,:);       

    % perform forward-backward kalman filtering
    KFTrack1 = kalmanFilterFwd(posData,kalmanOptions);
    KFTrack1 = kalmanFilterBkd(KFTrack1);        

    % flip the filters around to predict early spot positions in frames
    % prior to detection
    lastInd= min([5 length(frameVec)]);
    posDataTrunc = flipud(posData(1:find(framesFull==frameVec(lastInd),1),:));
    KFTrack2 = kalmanFilterFwd(posDataTrunc,kalmanOptions);
    KFTrack2 = kalmanFilterBkd(KFTrack2);
    addIndicesTo = 1:find(framesFull==frameVec(1),1)-1;
    addIndicesFrom = size(posDataTrunc,1):-1:size(posDataTrunc,1)-length(addIndicesTo)+1;
    
    % Add inferred position info to structure
    ParticlesTemp.framesFull = framesFull;
    ParticlesTemp.logL = KFTrack1.logL';
    ParticlesTemp.logLMean = nanmean(ParticlesTemp.logL);
    ParticlesTemp.xPosInf = KFTrack1.smoothedTrack(:,1);
    ParticlesTemp.yPosInf = KFTrack1.smoothedTrack(:,2);
    ParticlesTemp.zPosDetrendedInf = KFTrack1.smoothedTrack(:,3);    
    ParticlesTemp.xPosSEInf = sqrt(KFTrack1.smoothedTrackSE(:,1));
    ParticlesTemp.yPosSEInf = sqrt(KFTrack1.smoothedTrackSE(:,2));
    ParticlesTemp.zPosSEInf = sqrt(KFTrack1.smoothedTrackSE(:,3));
    
    % supplement with early point predictions
    ParticlesTemp.xPosInf(addIndicesTo) = KFTrack2.smoothedTrack(addIndicesFrom,1);
    ParticlesTemp.yPosInf(addIndicesTo) = KFTrack2.smoothedTrack(addIndicesFrom,2);
    ParticlesTemp.zPosDetrendedInf(addIndicesTo) = KFTrack2.smoothedTrack(addIndicesFrom,2);
    ParticlesTemp.zPosInf = ParticlesTemp.zPosDetrendedInf + trackingInfo.zPosStage(framesFull);
    ParticlesTemp.xPosSEInf(addIndicesTo) = sqrt(KFTrack2.smoothedTrackSE(addIndicesFrom,1));
    ParticlesTemp.yPosSEInf(addIndicesTo) = sqrt(KFTrack2.smoothedTrackSE(addIndicesFrom,2));
    ParticlesTemp.zPosSEInf(addIndicesTo) = sqrt(KFTrack2.smoothedTrackSE(addIndicesFrom,3));
    
    % add full prediction arrays as well
    ParticlesTemp.smoothedPredictions = KFTrack1.smoothedTrack;
    ParticlesTemp.smoothedPredictionsSE = sqrt(KFTrack1.smoothedTrackSE);
    
    ParticlesTemp.smoothedPredictions(addIndicesTo,:) = KFTrack2.smoothedTrack(addIndicesFrom,:);
    ParticlesTemp.smoothedPredictionsSE(addIndicesTo,:) = sqrt(KFTrack2.smoothedTrackSE(addIndicesFrom,:));
    
    % make filter for convenience
    ParticlesTemp.obsFrameFilter = ismember(framesFull,frameVec);