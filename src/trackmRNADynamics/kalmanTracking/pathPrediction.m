function ParticlesTemp = pathPrediction(ParticlesTemp, backwardTracks, trackingInfo, kalmanOptions)

  
    [frameVec, frameOrder] = sort(backwardTracks.Frame);
    
    % generate sorted position info
    framesFull = max([1,frameVec(1)-trackingInfo.nExtrapFrames]):min([trackingInfo.nFrames,frameVec(end)+trackingInfo.nExtrapFrames]);
    posData = NaN(length(framesFull), size(backwardTracks.MeasurementVec(:,1:3),2));
    posData(ismember(framesFull,frameVec),:) = backwardTracks.MeasurementVec(frameOrder,1:3);       

    % perform forward-backward kalman filtering
    KFTrack1 = kalmanFilterFwd(posData,kalmanOptions);
    KFTrack1 = kalmanFilterBkd(KFTrack1);        

    % flip the filters around to predict early spot positions in frames
    % prior to detection    
    posDataFlipped = flipud(posData);
    % 1:find(framesFull==frameVec(lastInd),1)
    KFTrack2 = kalmanFilterFwd(posDataFlipped,kalmanOptions);
    KFTrack2 = kalmanFilterBkd(KFTrack2);  
    
    % Add inferred position info to structure
    ParticlesTemp.framesFull = framesFull;
    ParticlesTemp.logL = nanmean([KFTrack1.logL flipud(KFTrack2.logL)],2); 
    ParticlesTemp.logLMean = nanmean(ParticlesTemp.logL);
    ParticlesTemp.logLArray = nanmean(cat(3,KFTrack1.logLArray, flipud(KFTrack2.logLArray)),3); 
    % add predictions
    ParticlesTemp.smoothedPredictions = nanmean(cat(3,KFTrack1.smoothedTrack ,flipud(KFTrack2.smoothedTrack)),3);
    ParticlesTemp.smoothedPredictionsSE = nanmean(cat(3,KFTrack1.smoothedTrackSE, flipud(KFTrack2.smoothedTrackSE)),3);
    ParticlesTemp.xPosInf = ParticlesTemp.smoothedPredictions(:,1);
    ParticlesTemp.yPosInf = ParticlesTemp.smoothedPredictions(:,2);
    ParticlesTemp.zPosDetrendedInf = ParticlesTemp.smoothedPredictions(:,3);
    ParticlesTemp.xPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,1);
    ParticlesTemp.yPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,2);
    ParticlesTemp.zPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,3);
    
    % supplement with early point predictions
    ParticlesTemp.zPosInf = ParticlesTemp.zPosDetrendedInf + trackingInfo.zPosStage(framesFull)';           
    
    % make filter for convenience
    ParticlesTemp.obsFrameFilter = ismember(framesFull,frameVec);
    
%     % lastly, calculate individual and overall particle likelihood scores
%     ParticlesTemp.logLDistance = (posData-ParticlesTemp.smoothedPredictions).^2./ParticlesTemp.smoothedPredictionsSE.^2;
